/*lessthans.hpp: This file is made to answer the following question: given two CTD's, ctd1 and ctd2, and
  one CGen, cg, what are all the CGens less than cg with respect to ctd1 and ctd2? When the index of cg
  is large, there are many possibilities to consider (and therefore potentially many less thans to
  generate), so generating a whole list at once can become too expensive. Instead, this file provides
  an LTGen object to generate the less thans one by one. Throughout the comments and function parameters,
  ctd1, ctd2, and cg will be exactly as above.

  Generating less thans is by far the most computationally intensive part of this program. Thankfully, we
  have found a few heuristics that can help reduce computation time. We allow customization of one of
  these heuristics (namely, memoization) by making it a policy for the LTGen object and using template parameters.
  This file also provides implementations of some basic memoization policies.*/
#pragma once

#include "data.hpp"

#include <vector>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/filesystem.hpp>

using namespace boost::filesystem;
using std::vector;

/************************************************************************************************************
  Memoizer Policy
 ***********************************************************************************************************
 A choice of memoizer policy will allow LTGen to record the results of less than searches, greatly speeding
 up processes that might require trying the same less than search many different times. Memoizer classes
 must expose the following methods.
 (1) void start(const CTD *ctd1, const CTD *ctd2, const CGen& cg)
 Prepares to memoize the less than search specified by the given parameters (see the declaration
 of LTGen below for more information on what these parameters are and how they specify a search). This
 is called by LTGen.start() and should initialize the memoizer.
 (2) bool hasNext()
 This function will only be called after a less than search is initialied with start(). It should return
 true if the memoizer has some or all of the results for the current search memoized and if not all
 of the memoized results have been returned by calls to nextLT() since the last call to start(). This is
 used by LTGen when generating the next search result to see if the result is already memoized.
 (3) bool hasAll()
 This function will only be called after a search is initialized with start(). It should return true if
 and only if we have all of the results for the current search memoized. LTGen will use this function in
 conjunction with hasNext() to differentiate between full, partial, and no memoization.
 (4) void memoizeLT(const CIP& path, int_fast16_t& pathind, int_fast16_t& numhs)
 This function will be called repeatedly after start() whenever hasNext() returns false (signalling that
 no results are memoized for the current less than search). It should memoize the single less than result
 given by the input parameters. (Here, path is the CIP of the less than CGen, pathind is the last index
 in the CIP that is in use to define the path of the CGen, and numhs is the number of edges labelled "h"
 on the path. See the documentation of the LTGen class for more information.) This
 function will be used by LTGen to memoize the results of the less than search as they are generated.
 (5) void nextLT(CIP& path, int_fast16_t pathind, int_fast16_t& numhs)
 This function will be called repeatedly after start() as long as hasLTs() returns true. It should place
 the path and number of hs of the next search result into the references passed to it. (By "next search
 result," we mean the entry in the memoized list of less than results that comes after the previous list
 entry given by a call to this function, or the first entry in the memoized list of results if no
 previous call to this function has been made since the last call to start()). This function will be
 used by LTGen to access each of the memoized results sequentially.
 (6) void end(bool gotall)
 This function is called by LTGen.end() to clean up after the search. The parameter gotall is true if all
 results for this search were generated. This parameter can be used by the memoizer to tell whether the
 results memoized in the search are full or partial.

TODO: Figure out how to handle I/O or other errors that occur during memoization.*/

/*The first implementation of a memoizer policy that we provide is the NoMemoizer policy, which does not
  memoize results. This is the default policy used by LTGen.*/
class NoMemoizer
{
    public:
        inline void start(const CTD *, const CTD *, const CGen&) {}
        inline bool hasNext() { return false; }
        inline bool hasAll() { return false; }
        inline void memoizeLT(const CIP&, int_fast16_t&, int_fast16_t&) {}
        inline void nextLT(CIP&, int_fast16_t&, int_fast16_t&) {}
        inline void end(bool) {}
};

/*The second implementation of a memoizer policy that we provide is the FSMemoizer policy, which memoizes
  the results of each less than search by storing them in their own file and uses a directory tree to easily
  find the file corresponding to the current search. This memoization schema has the advantage of leaving
  most of the complicated data management to the file system, which makes it portable and parallelizable.

  In terms of the particular structure of the directory tree that memoizes results, the base directory of
  the tree is given by FSMEM_BASE_DIR, which is defined below this comment. A less than search is
  specified precisely by the parameters passed to the start() function of the memoizer policy. We use
  CTD.toString() and CGen.toString() to turn these objects into strings. With this conversion, the results
  for the search will be memoized in FSMEM_BASE_DIR/ctd1/ctd2/cg.<ext>, where ctd1, ctd2, and cg denote
  the string representations of the paremeters of start() and where <ext> is "part" if the memo is
  an incomplete list of less than results and "full" if the memo has all the less than results.*/
const string FSMEM_BASE_DIR = "memos/FSMemoizer";
class FSMemoizer
{
    private:
        //In any single less than search, we will be using some memo file to read and/or write memoized results.
        //We store this file in an fstream so that the same variable can be used both when we read and when we
        //write.
        boost::filesystem::fstream file;
        //Because we may have to designate the memo file as full rather than partial at the end of the search,
        //which requires renaming the file, we also store the file name (without a .full or .part extension).
        path filename;
        //These two booleans indicate the state of memoization for this search. hassome is set to true
        //if we have partial or full results memoized for this search before beginning the search. hasall is set
        //to true if we have full results memoized for this search.
        bool hassome, hasall;
    public:
        void start(const CTD *ctd1, const CTD *ctd2, const CGen& cg);
        bool hasNext();
        bool hasAll();
        void memoizeLT(const CIP& path, int_fast16_t& pathind, int_fast16_t& numhs);
        void nextLT(CIP& path, int_fast16_t& pathind_out, int_fast16_t& numhs_out);
        void end(bool gotall);
        virtual ~FSMemoizer();
};

/************************************************************************************************************
  Generating Less Thans with LTGen
 ***********************************************************************************************************
 This class is our main tool for generating less thans. It generates them one-by-one, which allows us to use
 only one CIP allocation to represent all of the less than search results. One further simplification that
 we make is to consider less than results as being not CGens but rather pairs consisting of one CIP (the path
 of the CGen), and one integer (the number of edges of the path labelled "h"). In general, this can actually
 specify multiple CGens, as there may be multiple valid ways to label the edges. Storing all of these
 possible labellings separately would mean duplicating the same path over and over, which takes extra space
 and time to store (e.g. via a memoizer). Moreover, there are uses of LTGen which don't necessarily require
 choosing a valid labelling. (Perhaps the most important is when we use Theorem 1.19 to obstruct symplectic
 embeddings: here, we generate a possible decomposition by factorizing the input and finding less thans for
 each of the factorizations. We then label all of the CGens in the decomposition. This is more efficient
 because bullet 2 of Theorem 1.19 gives extra conditions that our labellings have to meet, which we can use
 to avoid ever considering certain label combinations.)

 The most usual pattern of use for LTGen is to iterate through all of the less thans of a given generator.
 The following loop shows the intended pattern of use of the methods of LTGen in order to accomplish this.

 LTGen lts;
 lts.start(ctd1, index, action, hcond); //Prepare the LTGen for the less than search.
 while (!lts.atEnd()) //As long as we haven't found the last one, keep generating less thans.
 {
      const CIP& path = lts.currPath(); //The first less than has already been generated. Get its path.
      int numhs = lts.currNumhs(); //Get the number of hs of the current less than.
      //Do stuff with path and numhs
      lts.next(); //Generate the next less than.
 }
 lts.end() //Let the generator clean up.
 //lts can now be used for another search -- just call start() again and repeat the above pattern.
 TODO: Should we make a class to wrap this pattern for the purpose of generating all less than CGens
 (possibly counting different labellings separately)?*/
template <class Memoizer = NoMemoizer>
class LTGen
{
    private:
        //This is the instance of the less than memoizer that we'll use.
        Memoizer memoizer;
        //These variables tell us what we're currently running a search on. More
        //precisely: if we want to find less thans for CGen cg with respect to
        //CTDs ctd1 and ctd2, then ctd1 here will just be ctd1, and action will
        //be the action of ctd2 on cg.
        const CTD* ctd1;
        double action;
        //Upper bound on the x-value and y-value of all less thans. These are
        //computed at the beginning of the search by CTD->LT_boundX() and CTD->LT_BoundY().
        int_fast16_t maxx, maxy;
        /*These variables hold information about the current state of the
          search. More precisely:
          - index is the "remaining index needed," i.e. the difference between the target index and
                  the index of the path we've built so far. When we find a valid less than, this
                  should be 0.
          - hcond is the "remainind h condition," i.e. the difference between x(cg) + y(cg) + m(cg) - 1
                  and the value of x+y - h/2 for the path we've built so far. When we find a valid less
                  than, this should be >= 0.
          - y is the y value of the current path we're trying.
          - remainingy is the y-value of the last vertex in the path we've built so far. When we find a
                       valid less than, this should be the length of the vertical edge (or 0 if none exists).
          - xsofar is the x-value of the last vertex in the path we've built so far. When we find a valid
                   less than, this will be the x-value of the entire path.
          - maxhs is the maximum number of h's we can put on the path we've built so far (in other words,
                  the total number of non-horizontal and non-vertical edges on the path).
          - topxs holds the maximum x-value for each edge. During the search, we first make an edge with
            the minimum possible x-value and then increase this x-value over and over; these values tell us
            when we can't increase the x-value anymore.
        */
        int index, hcond;
        int_fast16_t y, remainingy, xsofar, maxhs;
        vector<int_fast16_t> topxs;
        //These variables hold the current less than search result.
        CIP path;
        int_fast16_t pathind, numhs;
        //This variable tells us if we're done with the current search.
        bool done;
        //Helper methods to do the main tasks in next().
        inline void finishLT();
        inline bool addNewEdge();
        inline bool changeOrRemoveLastEdge();
        inline void finishBacktracking();
        //Helper methods used by the above ones to perform basic edge operations.
        inline bool tryAddingDiagonalEdge(int_fast16_t nexty);
        inline bool tryAddingHorizontalEdge(int_fast16_t length);
        inline bool tryIncreasingLastEdgeX();
        //More helper methods, used by the above ones for even more atomic operations.
        inline bool pathIsValid();
        inline void updateStatsForNewY();
        inline void updateStatsForNewDiagonalEdge();
        inline void updateStatsForNewHorizontalEdge();
        inline void updateStatsForRemovedEdge();
        inline int_fast16_t getDiagonalEdgeMinimumX(int_fast16_t nexty);
        inline int_fast16_t getDiagonalEdgeMaximumX(int_fast16_t nexty);
    public:
        //Accessors for the current less than search result.
        const CIP& currPath() const;
        int currLength() const;
        int currNumhs() const;
        //This function prepares the LTGen for the less than search specified by the given parameters. It
        //also generates the first less than, so make sure you use currPath() and currNumhs() to work with this
        //result before you call next()!
        void start(const CTD *ctd1, const CTD *ctd2, const CGen& cg);
        //This function returns true if and only if we have generated all of the results for this search.
        bool atEnd() const;
        //This function generates the next less than. If no next result exists, then done is set to true (so
        //that atEnd() will return true), and path and numhs are unchanged.
        void next();
        //This functions is called after the less than search is over to perform cleanup.
        void end();
};

/* Potential optimization of less-than algorithm
  1. Use arrays instead of vectors?
  2. Search edges in order of slope. Can possibly avoid GCD computations while doing this; edges
     can be generated using Farey sequences stuff (possibly with values of small Farey sequences hard-coded in
     or pre-generated on startup).
  3. If we don't avoid GCD computations with 2, then memoize GCD's (and/or pre-generate answers for low values).
  4. Compute actions incrementally as we search (to avoid recomputing contributions of a single edge.
     Most useful for LinearCTDs; may also speed up Ellipsoids and Polydisks, but unclear.)
  5. Pruning common subtrees of two subtrees (e.g. when two branches of search have identical nodes).
     Maybe this can be with a "transposition table," or as a dynamic programming algorithm?
*/
