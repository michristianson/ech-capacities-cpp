#define NODEBUG

#include "lessthans.hpp"

#include <string>
#include <cmath> // For ceil()
#include <iostream> // For cout (for debug print statements)

#include <boost/integer/common_factor.hpp>

using namespace std; // For cout (for debug print statements)

//We need to explicitly instantiate the versions of the template that we'll be using here.
//This solves the linking errors that come with splitting the declaration and definition of
//a template class into two files.
template class LTGen<NoMemoizer>;
template class LTGen<FSMemoizer>;

void FSMemoizer::start(const CTD* ctd1, const CTD* ctd2, const CGen& cg)
{
    string name = FSMEM_BASE_DIR + "/" + ctd1->toString()
                                 + "/" + ctd2->toString()
                                 + "/" + cg.toString();
    filename = current_path() / name;
    if (exists(filename.replace_extension("full")))
    {
        //We have full results memoized.
        hassome = true;
        hasall = true;
        file.open(filename, ios::in);
    }
    else if (exists(filename.replace_extension("part")))
    {
        //We have only a partial list of results memoized.
        file.open(filename);
        hassome = true;
        hasall = false;
    }
    else
    {
        //We have no results memoized, so we prepare to start memoizing results.
        hassome = false;
        hasall = false;
        //Create the directories that the memo file will live in.
        create_directories(filename.parent_path());
        file.open(filename, ios::out);
    }
}
bool FSMemoizer::hasNext()
{
    //If we didn't have any results memoized for this search, then there can't be a next result to read.
    if (!hassome)
        return false;
    //To see if we're at the end of a memo, we peek at the next character of a file and
    //see if this sets the eof flag (which will indicate that we're at the end of the file).
    file.peek();
    return !file.eof();
}
bool FSMemoizer::hasAll() { return hasall; }
void FSMemoizer::memoizeLT(const CIP& path, int_fast16_t& pathind, int_fast16_t& numhs)
{
    /*We will essentially serialize the parameters into the memo file. We put each less than
      on its own line; this doesn't actually help us read them from the memo any more easily,
      but it makes our memo files easier for human eyes to parse. Our serialization schema will
      serialize a CIP of length n with edges (x_1,y_1), (x_2,y_2), ..., (x_n,y_n) as
n:x1,y1;x2,y2;...;xn,yn;numhs*/
    file << (pathind + 1) << ":";
    for (int_fast16_t i = 0; i <= pathind; i++)
        file << path[i].first << "," << path[i].second << ";";
    file << numhs << "\n";
}
void FSMemoizer::nextLT(CIP& path, int_fast16_t& pathind_out, int_fast16_t& numhs_out)
{
    //According to the serialization schema described in memoizeLT(), the first thing in the
    //serialization we need to read is the number of edges of the CIP of the less thans.
    file >> pathind_out;
    //Skip the colon after the file size.
    file.ignore(1);
    //Resize the path so that it has the right number of elements.
    path.resize(pathind_out);
    //Now, iterate over all the edges, adding them to the CIP as we go.
    int_fast16_t x, y;
    for (size_t i = 0; i < pathind_out; i++)
    {
        //Get the x-coordinate of the edge.
        file >> x;
        //Skip the comma separating the coordinates.
        file.ignore(1);
        //Get the y-coordinate of the edges.
        file >> y;
        //Store the edge in the CIP.
        path[i] = make_pair(x, y);
        //Skip the semicolon separating edges. Note that there is also a semicolon after the last edge
        //but before the number of hs, and we need to skip that one too, so it is good to do this even
        //in the last iteration of the loop.
        file.ignore(1);
    }
    //The last thing to do is to read out the number of hs.
    file >> numhs_out;
    //Skip the newline at the end of this less than result.
    file.ignore(1);
}
void FSMemoizer::end(bool gotall)
{
    //Since we're done, we can close the memo file.
    file.close();
    //If we didn't already have all the results memoized before this search, then the results
    //are currently stored in a file with a .part extension. If we generated all the results
    //with this search, then we should mark this memo as full rather than partial.
    if (!hasall && gotall)
    {
        path newname = filename;
        newname.replace_extension("full");
        rename(filename, newname);
    }
}
FSMemoizer::~FSMemoizer() { file.close(); }


//LTGen<Memoizer> accessor methods for the current less than.
template <class Memoizer>
const CIP& LTGen<Memoizer>::currPath() const { return path; }
template <class Memoizer>
int LTGen<Memoizer>::currLength() const { return pathind + 1; }
template <class Memoizer>
int LTGen<Memoizer>::currNumhs() const { return numhs; }


//LTGen<Memoizer> methods to perform less than search.
template <class Memoizer>
void LTGen<Memoizer>::start(const CTD *ctd1, const CTD *ctd2, const CGen& cg)
{
#ifdef DEBUG
    cout << "Initializing less than search." << endl;
#endif
    this->ctd1 = ctd1;
    //We don't actually need to store ctd2 or cg; we just need to use them
    //to initialize some statistics relevant to the definition of "less than."
    index = cg.I();
    action = ctd2->actionOf(cg);
    hcond = cg.x() + cg.y() + cg.m() - 1;
    maxx = ctd1->LT_boundX(index, action, hcond);
    maxy = ctd1->LT_boundY(index, action, hcond);
#ifdef DEBUG
    cout << "Computed initial stats: index = " << index << ", action = " << action
         << ", hcond = " << hcond << endl;
    cout << "Upper bounds for search: maxx = " << maxx << ", maxy = " << maxy << endl;
#endif
    if (maxy < 0)
    {
        //If the maximum value of y is negative, then no possible y-value
        //exists, so there are no less thans.
        done = true;
        return;
    }

    //Initialize the memoizer.
    memoizer.start(ctd1, ctd2, cg);
    //Initialize the search variables.
    done = false;
    y = 0;
    remainingy = 0;
    xsofar = 0;
    maxhs = 0;
    /*Initialize the values that will contain the search results. Since each edge
      (except possibly one horizontal one) contributes at least 1 to the y-value of the CGen,
      we can only have at most maxy+1 edges. We resize path and topxs to this
      maximum length, so that we don't have to expand the size of these vectors on-the-fly.*/
    path.resize(maxy + 1);
    topxs.resize(maxy + 1);
    pathind = -1;
    numhs = 0; 

    /* Now we need to start the search itself. The first path Gen to try is a
     * completely horizontal one, so we go ahead and make this. (next() would
     * generate it automatically, but it would start with a horizontal path of
     * length 1, then extend it to length 2, and so on. We can skip this
     * extension process by just making the a completely horizontal path that
     * is "long enough" to begin with.)
     *
     * For a horizontal CGen of length x, we will have 
     *
     * I = 2(L-1) = 2x.

     * This means that x (hence the entire CGen) is completely determined by the
     * index condition of being a less than: we want x = index/2. In case index
     * is odd, we'll just take the floor. This makes the computations easy,
     * because floor(index / 2) = index >> 1.*/
#ifdef DEBUG
    cout << "Search initialized. Trying completely horizontal path." << endl;
#endif
    if (tryAddingHorizontalEdge(index >> 1) && pathIsValid())
        finishLT();
    else
        next();
}

template<class Memoizer>
void LTGen<Memoizer>::next()
{
#ifdef DEBUG
    cout << "Searching for next valid less than." << endl;
#endif
    // If the search is done, there's nothing to do.
    if (done)
    {
#ifdef DEBUG
        cout << "Search is done. Nothing else to do." << endl;
#endif
        return;
    }

    //If we have a result memoized, use that.
    if (memoizer.hasNext())
    {
#ifdef DEBUG
        cout << "Using memoizer to get next less than." << endl;
#endif
        memoizer.nextLT(path, pathind, numhs);
        return;
    }
    //If we don't have a result memoized, but the memorizer has all the results,
    //then we've already returned all results (using the memoizer), so we're done.
    else if (memoizer.hasAll())
    {
#ifdef DEBUG
        cout << "Memoizer has returned all search results, so the search is done." << endl;
#endif
        done = true;
        return;
    }

    /*The search assumes there are no vertical edges in the path; so, if we
      returned a path with a vertical edge on the last call to next(), we need
      to remove it. Fortunately, vertical edges change none of our statistics,
      so all we have to do is remove the vertical edge from the current path.*/
    if (pathind >= 0 && path[pathind].first == 0)
    {
#ifdef DEBUG
        cout << "Vertical edge removed from path for search." << endl;
#endif
        pathind--;
    }

    bool justbacktracked = false;
    //Here we unroll the first iteration of the while loop below, with two small
    //adjustments:
    //(1) We don't check if we've found a valid path, since we need to actually
    //    change the path at least once before we've found the next valid path.
    //(2) We already know that justbacktracked = false the first time around,
    //    which makes things a bit simpler.
    if (!addNewEdge())
    {
        if (!changeOrRemoveLastEdge())
        {
            finishBacktracking();
            if (done)
                return;
            if (pathind >= 0)
                justbacktracked = true;
        }
    }
    /*The while (true) loop below essentially serves to "flatten" a simple
      recursive algorithm to search for all possible paths. The recursive algorithm
      would look something like this (note that we use "yield" and "yield from"
      like the Python keywords for generators, which mean "return this result,
      and then pick up again from here the next time we call next()):

      next()
      {
          for (y = 0; y <= maxy; y++) { yield from recursiveSearch(); }
      }
      recursiveSearch()
      {
          if (pathIsValid()) { yield path; }
          while (hasNextEdgeToTry())
          {
              changeToNextEdge();
              yield from recursiveSearch();
          }
          removeEdgeFromPath();
      }

      In our "flattened" version below:
      - the if(pathIsValid()) part remains largely the same;
      - addnewEdge() is the first iteration of the while loop in recursiveSearch()
        (including the check of whether there are any edges to try in the first place);
      - changeOrRemoveLastEdge() takes care of all subsequent iterations of the while
        loop as well as removing the edge at the end of recursiveSearch();
      - finishBacktracking() serves to "flatten out" the for loop on y in our
        above version of next().

      The bool justbacktracked in our "flattened" version keeps track of whether
      we are starting a new recursion (in which case we should check for a valid
      path, and if not, then do the first iteration of the while loop)) or
      whether we just ended a recursive call and backtracked to the caller
      (in which case we need to change to the next edge to continue the while
      loop in the caller, or remove the edge if we've finished the while loop).
      */
    while(true)
    {
#ifdef DEBUG
        cout << "Entering while (true) loop. Path currently has length " << (pathind+1)
             << ", y-value " << y << ", and x-value " << xsofar
             << ", and justbacktracked = " << (justbacktracked ? "true." : "false.") << endl;
        if (pathind >= 0)
        {
            cout << "Current path: { ";
            for (int i = 0; i < pathind; i++)
                cout << "(" << path[i].first << "," << path[i].second << "), ";
            cout << "(" << path[pathind].first << "," << path[pathind].second << ") }" << endl;
        }
#endif
        //If we just backtracked, we've already tried adding a new edge here.
        //Otherwise, it's time to try adding a new edge.
        if (!justbacktracked)
        {
            /*If we've found a valid path, we're done. We only check this when
              we haven't just backtracked, because otherwise we could backtrack
              to a valid path that we've already had. (More specitifically: if
              we return a path ending in a vertical edge, then call next() again
              and try a bunch of non-vertical edges instead; when we backtrack
              past those edges, the search will see a "valid path" again,
              because we're back to the path we already returned minus the
              vertical edge.)*/
            if (pathIsValid())
            {
                finishLT();
                return;
            }
            if (addNewEdge())
            {
                //After this, we will have:
                //justbacktracked = false;
                continue;
            }
        }
        //We only get here if we just backtracked or we tried to add a new
        //edge in above if, but we couldn't. Either way, there are no more
        //new edges to try adding; instead, try changing the current edge if
        //possible, and remove it otherwise.
        if (changeOrRemoveLastEdge())
        {
#ifdef DEBUG
            cout << "Edge changed successfully; setting justbacktracked = false." << endl;
#endif
            justbacktracked = false;
        }
        else
        {
            finishBacktracking();
            if (done)
                return;
            //If we just backtracked past all the edges, we want to try adding
            //a new edge, not changing the last edge. In a sense, backtracking
            //past the last edge means we didn't just backtrack: instead, we
            //just increased the y-value.
            if (pathind < 0)
            {
#ifdef DEBUG
                cout << "Just backtracked past first edge, so setting justbacktracked = false anyway." << endl;
#endif
                justbacktracked = false;
            }
            else
            {
#ifdef DEBUG
                cout << "Just backtracked, so setting justbacktracked = true." << endl;
#endif
                justbacktracked = true;
            }
        }
    }
}



//Helper methods to perform the main tasks in LTGen<Memoizer>::next().
template<class Memoizer>
inline void LTGen<Memoizer>::finishLT()
{
#ifdef DEBUG
    cout << "Found valid less than; finishing it up now." << endl;
#endif
    if (remainingy != 0)
    {
        //Add a vertical edge to get us down to the x-axis.
        pathind++;
        path[pathind].first = 0;
        path[pathind].second = -remainingy;
    }
    //Set the number of hs such that we get the correct index.
    numhs = -index;
    memoizer.memoizeLT(path, pathind, numhs);
}
template<class Memoizer>
inline bool LTGen<Memoizer>::addNewEdge()
{
#ifdef DEBUG
    cout << "Trying to add a new edge at pathind = " << (pathind+1) << endl;
#endif

    for (int_fast16_t nexty = -remainingy; nexty < 0; nexty++)
    {
        if (tryAddingDiagonalEdge(nexty))
            return true;
    }
    /*The last edge y-value to try is 0, i.e. a horizontal edge. This is
      only possible when there are no edges already in the path; in this
      case, we try adding in a horizontal edge of the minimum length, which
      is 1.*/
    if (pathind < 0 && tryAddingHorizontalEdge(1))
        return true;
    return false;
}
template<class Memoizer>
inline bool LTGen<Memoizer>::changeOrRemoveLastEdge()
{
#ifdef DEBUG
    cout << "Trying to change or remove edge at pathind = " << pathind << endl;
#endif

    //If there is no last edge, then there's nothing to do.
    if (pathind < 0)
        return false;

    //First, try increasing the x-value of the last edge.
    if (tryIncreasingLastEdgeX())
    {
        return true;
    }

    int_fast16_t nexty = path[pathind].second + 1;
    updateStatsForRemovedEdge();
    pathind--;
    /*If we can't increase the x-value of the edge anymore, then we need to
      increase the y-value and reset the x-value. To do this, we remove the
      last edge (which we did above if it exists) and then behave as if we're
      trying to add a new edge back on with a smaller y-value.
      
      Note: if the last edge was horizontal, then nexty == 1 here, and there
      are no more edges to try, so we avoid that case here. */
    if (nexty < 1)
    {
        //Try all smaller y-values until we find one that works.
        //TODO: It might be the case that, once one y-value doesn't work, all
        //smaller y-values don't work. If so, then we don't need to loop here;
        //instead, we can just try the very first value of nexty.
        while (nexty < 0)
        {
            if (tryAddingDiagonalEdge(nexty))
                return true;
            else
                nexty++;
        }
        /*The last edge y-value to try is 0, i.e. a horizontal edge. This is
          only possible when there are no edges already in the path; in this
          case, we try adding in a horizontal edge of the minimum length, which
          is 1.*/
        if (pathind < 0 && tryAddingHorizontalEdge(1))
            return true;
    }
#ifdef DEBUG
    cout << "Could not change last edge, so it was removed. " << endl;
#endif
    return false;
}
template<class Memoizer>
inline void LTGen<Memoizer>::finishBacktracking()
{
    //Note that changeOrRemoveLastEdge() already removed the last edge if there
    //was no new edge to change it to. So, all we have to do here is bump up the
    //y-value of the path if we've backtracked past all the edges in the path.
    if (pathind < 0)
    {
        y++;
#ifdef DEBUG
        cout << "---------------------------------" << endl;
        cout << "y-value of path increased to " << y << endl;
        cout << "---------------------------------" << endl;
#endif
        /*Note: it's okay for index to be 0 AFTER updating stats for the new
          y-value: this just means a single vertical edge is the only possible
          path. But if index == 0 right now, then index < 0 after we update it
          with the new y-value, and that will make it impossible to get anything
          with the correct index (adding any edge adds at least 2 to the index but
          only provides at most 1 more edge labelled 'h' to bring the index back down).*/
        if (y > maxy || index == 0)
        {
            done  = true;
#ifdef DEBUG
            cout << "y-value too large now; search must be done." << endl;
#endif
        }
        else
        {
            updateStatsForNewY();
        }
    }
}




//Helper methods for LTGen<Memoizer>::next() to perform basic operations on
//edges in the path.
template<class Memoizer>
inline bool LTGen<Memoizer>::tryAddingDiagonalEdge(int_fast16_t nexty)
{
#ifdef DEBUG
    cout << "Trying to add diagonal edge with y-value " << nexty << endl;
#endif
    //If remainingy == 0 or the index is already too low, then there's
    //no more room to add a new diagonal edge.
    if (remainingy == 0 || index < -maxhs)
    {
#ifdef DEBUG
        cout << "No more room to add diagonal edge" << endl;
#endif
        return false;
    }

    //Find the x-value of the new edge, make sure it's valid.
    path[pathind + 1].first = getDiagonalEdgeMinimumX(nexty);
    topxs[pathind + 1] = getDiagonalEdgeMaximumX(nexty);
    if (path[pathind + 1].first > topxs[pathind + 1])
    {
#ifdef DEBUG
        cout << "Failed to add edge because x-values are too large" << endl;
#endif
        return false;
    }

    //At this point, we know we can add a new edge successfully. We just need to
    //finish storing the edge in path and then update all of our statistics.
    pathind++;
    path[pathind].second = nexty;
    updateStatsForNewDiagonalEdge();
#ifdef DEBUG
    cout << "Successfully added edge with x-value " << path[pathind].first << endl;
#endif
    return true;
}
template<class Memoizer>
inline bool LTGen<Memoizer>::tryAddingHorizontalEdge(int_fast16_t length)
{
    //Since horizontal edges are the first edge in the path, the only way to
    //fail to add one is if it's too long for maxx.
    //TODO: Should we also check that the index isn't too large, like in
    //tryAddingDiagonalEdge() and tryIncreasingLastEdgeX()? This might
    //avoid some backtracking on horizontal edges that are clearly too long.
    if (length > maxx)
    {
#ifdef DEBUG
        cout << "Failed to add horizontal edge because length > maxx." << endl;
#endif
        return false;
    }
    pathind++;
    path[pathind].first = length;
    path[pathind].second = 0;
    topxs[pathind] = maxx;
    updateStatsForNewHorizontalEdge();
#ifdef DEBUG
    cout << "Horizontal edge of length " << length << " added successfully." << endl;
#endif
    return true;
}
template<class Memoizer>
inline bool LTGen<Memoizer>::tryIncreasingLastEdgeX()
{
#ifdef DEBUG
    cout << "Trying to increase x-value of last edge from " << path[pathind].first
         << ". Bound from topxs is " << topxs[pathind] << "." << endl;
#endif
    //We can only increase the x-value of the last edge if it isn't already at
    //its upper bound.
    if (path[pathind].first >= topxs[pathind])
    {
#ifdef DEBUG
        cout << "Could not increase x-value because of topxs bound." << endl;
#endif
        return false;
    }

    /* This equation for how the index will change is just newindex = index - change
     * with
     *
     * change (index update for (first+1,second)) - (index update for (first, second)).
     *
     * Here, index update is just the equation used in updateStatsForNewEdge() to
     * update index. There's just one catch: since the edge has already been
     * added, remainingy has changed from its value in updateStatsForNewEdge().
     * The value there is the value it has here, minus path[pathind].second.
     */
    int newindex = index - 2 * (remainingy - path[pathind].second) - path[pathind].second - 1
        - boost::integer::gcd(path[pathind].first + 1, path[pathind].second)
        + boost::integer::gcd(path[pathind].first, path[pathind].second);
    if (newindex < -maxhs)
    {
        /*Increasing the x makes the index of the less than too large. Note that
          this restriction cannot be bundled into topxs, because we don't have a
          good way to determine ahead of time what will be the maximum x-value
          before the index gets too large. (The main issue is the appearance
          of a gcd in the index contribution of the edge, which makes it impossible
          to solve for the x-value of the edge explicitly in the equation
          newindex = -maxhs.)*/
#ifdef DEBUG
        cout << "Could not increase x-value because it would make the index too large." << endl;
#endif
        return false;
    }

    //Now we know that we can increase the x-value of the edge. We just need to
    //go ahead and update both the path and our stats to reflect this change.
    path[pathind].first++;
    index = newindex;
    hcond -= 1;
    xsofar += 1;
#ifdef DEBUG
    cout << "Successfully increased x-value of edge." << endl;
#endif
    return true;
}



// More helper methods for LTGen<Memoizer>::next().
template <class Memoizer>
inline bool LTGen<Memoizer>::pathIsValid()
{
#ifdef DEBUG
    cout << "Checking if path is valid: index = " << index << ", hcond = " << hcond
         << ", action = " << action << ", and action of proposed less than is "
         << ctd1->LT_actionOf(path, pathind, remainingy, xsofar, y) << endl;
#endif
    return index <= 0 && index >= -maxhs &&
        hcond <= index / 2 &&
        ctd1->LT_actionOf(path, pathind, remainingy, xsofar, y) <= action;
}
template<class Memoizer>
inline void LTGen<Memoizer>::updateStatsForNewY()
{
    index -= 2;
    hcond -= 1;
    remainingy++;
}
template<class Memoizer>
inline void LTGen<Memoizer>::updateStatsForNewDiagonalEdge()
{
    //Method to update stats on the current path (specifically, index, hcond,
    //xsofar, remainingy, and maxhs) when a new edge is added. We assume that
    //the new edge is already at path[pathind].

    /* We use Pick's Theorem here to compute the index contribution of the new
     * edge. This theorem tells us that the number of lattice points beneath
     * the edge is
     *
     * L = A + b/2 + 1, 
     *
     * where A is the area underneath the edge and b is the number of lattice
     * points on the boundary of the trapezoid formed by the edge, the x-axis,
     * and two vertical line segments from each end of the edge to the x-axis.
     * (See the comments in CGen::L() in data.cpp for more details.) Note that
     * the vertices of this trapezoid are (xsofar, 0), (xsofar + edge.first, 0),
     * (xsofar, remainingy), and (xsofar + edge.first, remainingy +
     * edge.second). Using this, we get
     *
     * A = 1/2 * (2 * remainingy + edge.second) * edge.first
     *
     * and
     *
     * b = edge.first + remainingy + remainingy + edge.second + gcd(edge.first, edge.second).
     *
     * (Note that gcd(edge.first, edge.second) + 1 is the number of lattice
     * points on the edge itself. The + 1 disappears in the above equation
     * because of "double-counting" corners.)
     *
     * Now, all of the lattice points on the left vertical line of the trapezoid
     * were already underneath the path (they all occur beneath the last vertex
     * in the path prior to adding this edge). There are (remainingy + 1) of
     * these. So, we subtract those off from L, and then we use I = 2L + (constants)
     * to find that this new edge adds the following value to the index of our path:
     *
     * 2*(L - remainingy - 1) = 2*A + b - 2*remainingy.
     *
     * The righthand side of this equation simplifies to
     * 
     * (2 * remainingy + edge.second + 1) * edge.first + edge.second + gcd(edge.first, edge.second)
     *
     * To update index, then, we just subtract this equation from it.
     */
    index = index - (2 * remainingy + path[pathind].second + 1) * path[pathind].first
                  - path[pathind].second
                  - boost::integer::gcd(path[pathind].first, path[pathind].second);
    // We update hcond by subtracting off changes to x and y. Adding a new edge
    // doesn't change y, but it does increase xx.
    hcond -= path[pathind].first;

    xsofar += path[pathind].first;
    //Note: path[pathind].second <= 0, so the below line makes remainingy smaller.
    remainingy += path[pathind].second;
    //The new edge can be labelled h provided it isn't horizontal (or vertical,
    //but we never add in vertical edges until we've got a less than to return,
    //so we should never be dealing with a vertical edge here).
    if (path[pathind].second != 0)
        maxhs++;
}
template<class Memoizer>
inline void LTGen<Memoizer>::updateStatsForNewHorizontalEdge()
{
    /*This method is the same as updateStatsForNewDiagonalEdge(), but specialized
      to the case where the edge is horizontal. Everything in
      updateStatsForNewDiagonalEdge() would be correct for a horizontal edge,
      but in the horizontal case we can simplify all of the equations because
      path[pathind].second == 0. */
    index = index - 2 * (remainingy + 1) * path[pathind].first;
    hcond -= path[pathind].first;
    xsofar += path[pathind].first;
    //Note: we don't update remainingy or maxhs, because they don't change when
    //we add a horizotal edge.
}
template<class Memoizer>
inline void LTGen<Memoizer>::updateStatsForRemovedEdge()
{
    /* Method to update stats after an edge is removed from the path.
     * Note: We assume in this function that the removed edge is still at path[pathind].
     * This means that, when removing the edge, one should call this function first,
     * then do pathind--.*/

    /* We essentially undo everything done in updateStatsForNewEdge(). It is important
     * that we reverse the order of the index update and the remainingy update when we
     * do this: the formula for updating updateStatsForNewEdge() is based on what
     * remainingy was before that update, so we need to reset it to this value first,
     * and then we can undo the index update.*/
    hcond += path[pathind].first;

    xsofar -= path[pathind].first;
    remainingy -= path[pathind].second;
    if (path[pathind].second != 0)
        maxhs--;

    index = index + (2 * remainingy + path[pathind].second + 1) * path[pathind].first
                  + path[pathind].second
                  + boost::integer::gcd(path[pathind].first, path[pathind].second);
}
template<class Memoizer>
inline int_fast16_t LTGen<Memoizer>::getDiagonalEdgeMinimumX(int_fast16_t nexty)
{
    /* Suppose we want to add edge (nextx, nexty) to the path. We want to make
     * sure that nextx is big enough (i.e. the path has small enough slope) that
     * we can still extend out path out to a path that has large enough index.
     * Since every edge after this one must have a steeper slope (to make a
     * convex path), an upper bound on the index is given by extending the edge
     * (nextx, nexty) all the way down to the x-axis in a straight line. This
     * means instead adding the edge (c*nextx, c*nexty), where c = -remainingy/nexty.
     *
     * We can compute the index contribution of this edge just like we did
     * in updateStatsForNewEdge(). The only difference is that c might not be an
     * integer, which complicates the lattice point calculations just slightly.
     * For the right triangle bounded by this new edge, the x-axis, and a vertical line
     * on the left, we get
     *
     * A = 1/2 * c^2 * nextx * (-nexty),
     * b = floor(c * nextx) + remainingy + gcd(nextx, -nexty) * floor(c).
     *
     * (The only difference between this and updateStatsForNewEdge() is the
     * floor functions, which come from the fact that c isn't necessarily an integer --
     * plus, of course, replacing edge.first and edge.second with c*nextx and c*nexty = -remainingy.)
     *
     * Using Pick's theorem and I = 2*L + (constants)  and subtracting off the
     * (remainingy + 1) points on the left edge of the right triangle, we see
     * that the index contribution of this edge is
     *
     * c^2 * nextx * (-nexty) + floor(c * nextx) + gcd(nextx, -nexty) * floor(c) - remainingy.
     *
     * Since we've extended our edge all the way to the x-axis, we can't add any
     * more edges, so we need to make sure the above value is at least as much
     * as the rest of the index contribution that we need. In other words, we
     * need
     *
     * index <= c^2 * nextx * (-nexty) + floor(c * nextx) + gcd(nextx, -nexty) * floor(c) - remainingy.
     *
     * Removing the floor functions makes this inequality still true, which gives
     *
     * index + remainingy <= nextx * (c^2 * (-nexty) + c) + c * gcd(nextx, -nexty).
     *
     * Solving for nextx, we get
     *
     * nextx >= [index + remainingy - c * gcd(nextx, -nexty)] / [(c^2 * (-nexty) + c)].
     *
     * Since, gcd(nextx, -nexty) <= (-nexty), we may replace the gcd with nexty
     * to get a valid inequality (making the negative gcd term larger will make the
     * righthand side smaller). Using the fact that c * (-nexty) = remainingy, we
     * then get
     *
     * nextx >= index / [c (remainingy + 1)]
     *       = index * (-nexty) / [remainingy * (remainingy + 1)]
     */
#ifdef DEBUG
    cout << "Computing lower bound on x for diagonal edge." << endl;
#endif

    /*Note: Here (and on hbound below), we use round() when it would
      mathematically make sense to use ceil() (since we are getting a lower
      bound on an integer). I'm not sure if this is necessary, but I think that
      round() will be more robust to floating-point error, and it is better to
      return a slightly weaker bound than to return an incorrect bound and fail
      to try some edges that we should try. (I did the same thing in
      LT_boundXUsingAction() and LT_boundYUsingAction() for Polydisk and
      Ellipsoid in data.cpp, because I actually found a situation where using
      floor() instead of round() gives the wrong answer there.) */
    int_fast16_t indexbound = round( index * (-nexty) / (remainingy * (remainingy + 1)) );
    /* We can also get a bound on x from the "h condition" of a less than. We
     * need the final x-value x of the path to satisfy hcond - x <= 0. As above,
     * if we plan to add edge (nextx, nexty), then because later edges have
     * steeper slopes, the largest x-value we can get is if we extend that edge
     * all the way to the y-axis, in which case we will add c * nextx to the
     * x-value so far, where c = -remainingy/nexty as above. This gives us
     *
     * c * nextx >= x >= hcond,
     *
     * hence
     *
     * nextx >= hcond / c = hcond * (-nexty) / remainingy.
     */
    int_fast16_t hbound = round( hcond * (-nexty) / remainingy );
    // In case the above bounds are not helpful, we make sure that the smallest
    // x-value we consider is 1 (vertical edges, i.e. those with x-value 0, are
    // only added at the end of the path if we've found a valid path).
#ifdef DEBUG
    cout << "Bounds found: indexbound = " << indexbound << ", hbound = " << hbound << endl;
#endif
    return max<int_fast16_t>(1, max<int_fast16_t>(indexbound, hbound));
}
template<class Memoizer>
inline int_fast16_t LTGen<Memoizer>::getDiagonalEdgeMaximumX(int_fast16_t nexty)
{
    /*We have two upper bounds on the x-value of an edge: the first is to ensure
      that the slope of the edge is steeper (i.e. more negative) than the
      previous edge, and the second is that it doesn't send the x-value of the
      CGen over the maximum possible x-value. The latter is simply maxx -
      xsofar; the former is the following inequality:

      nexty / nextx <= path[pathind].second / path[pathind].first.

      Note that if there is no previous edge, or if the previous edge is horizontal,
      then this bound based on slopes does not apply.

      Note: we assume that pathind has NOT been increased to account for the new
      edge yet. In other words, callers that aim to create a newedge should use pathind++
      AFTER calling this method. */
    if (pathind < 0 || path[pathind].second == 0)
    {
#ifdef DEBUG
        cout << "Computing upper bound on x for first non-horizontal edge in path. Bound is: " << (maxx-xsofar) << endl;
#endif
        return maxx - xsofar;
    }
    else
    {
#ifdef DEBUG
        cout << "Computing upper bound on x for new diagonal edge with nexty = " << nexty
             << ". Previous edge: (" << path[pathind].first << "," << path[pathind].second << ")." << endl;
#endif
        /*Note: For any number x, ceil(x) - 1 is the largest integer strictly
          less than x (floor(x), by contrast, would just give x if x is an
          integer, not something strictly smaller than x). This is important
          because our slopes should be strictly decreasing.

          TODO: Should we be worried about the possibility of floating-point error
          in the equation below? Thanks to ceil(), we could end up getting the bound 
          to be 1 less than it should be, which could mean we fail to try 1 value of
          x that we should have tried. We switched to using round() in
          getDiagonalEdgeMinimumX() just in case, but I don't think that would
          work here: the slope bound has no margin for error on either side, or
          else we risk missing an edge we should try or trying an edge with the
          same slope as the previous one. (By contrast, our lower bounds on x
          are still okay if they're a little too low -- we'll just try some
          extra x-values that won't work.)

          If this is a problem, I'm not quite sure what we could do to fix it.
          One failsafe would be to explicitly check if (upper bound on x)+1 is
          actually a valid possibility (either here or in tryIncreasingLastEdgeX()),
          but given that I haven't seen any problems here, this feels a little unnecessary.
        */
        int_fast16_t slopebound = ceil((double)path[pathind].first / (double)path[pathind].second * nexty) - 1;
#ifdef DEBUG
        cout << "Bounds found: slopebound = " << slopebound << ", maxbound = " << (maxx-xsofar) << endl;
#endif
        return min(maxx - xsofar, slopebound);
    }
}



template <class Memoizer>
bool LTGen<Memoizer>::atEnd() const
{ 
    return done;
}
template<class Memoizer>
void LTGen<Memoizer>::end()
{
    memoizer.end(done);
}
