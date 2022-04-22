#include "data.hpp"

#include <sstream>
#include <limits>

/* Note: After looking at the unit tests in boost's common_factor_tests.cpp,
 * I am fairly certain that boost::integer::gcd() has the following properties
 * for all integers m and n:
 *
 * (1) gcd(0,n) = gcd(n,0) = n.
 * (2) gcd(m,n) = gcd(-m,n) = gcd(m,-n) = gcd(m,n).
 *
 * These properties may not be entirely standard, but they do seem to be the
 * most logical default behavior in some way. They are also exactly the behavior
 * we want (both here and in lessthans.cpp) when we find gcd's of edges.
 */
#include <boost/integer/common_factor.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/bind/bind.hpp>


//Vertices and Edges
Edge operator+(const Edge& v1, const Edge& v2)
{
    return make_pair(v1.first + v2.first, v1.second + v2.second);
}
Edge operator-(const Edge& v1, const Edge& v2)
{
    return make_pair(v1.first - v2.first, v1.second - v2.second);
}
/*Check for IEC 559 (a.k.a. IEEE 754) arithmetic. This should be present on
  most modern machines; if it is, then division by 0 is allowed and will return
  INFINITY, -INFINITY, or NaN as appropriate. If not, then we attempt a quick
  workaround for calculating slopes. This workaround has not been tested
  extensively, but preliminary tests suggest that it probably works. */
#if __STDC_IEC_559__
    inline double edgeSlope(const Edge& e)
    {
        return (double)e.second / (double)e.first;
    }
    inline double vertexSlope(const Vertex& v1, const Vertex& v2)
    {
        return (double)(v2.second - v1.second) / (double)(v2.first - v1.first);
    }
#else
    #define VERTICAL_SLOPE -numeric_limits<double>::max()
    #warning "IEEE 754 arithmetic not detected, so division by 0 is not allowed. Attempting to fix this by viewing vertical lines as having slope -(max double value), but this approach has not be tested extensively and may break. Using a machine with IEEE754 arithmetic is recommended."
    inline double edgeSlope(const Edge& e)
    {
        if (e.first == 0)
        {
            if (e.second == 0)
                throw invalid_argument("Taking slope of edge (0,0)");
            return VERTICAL_SLOPE;
        }
        return (double)e.second / (double)e.first;
    }
    inline double vertexSlope(const Vertex& v1, const Vertex& v2)
    {
        if (v1.first == v2.first)
        {
            if (v1.second == v2.second)
                throw invalid_argument("Taking slope between two equal vertices.");
            //TODO: Should we be worried about which vertex comes first here to
            //determine if the caller actually wants -VERTICAL_SLOPE? I don't
            //think this should ever come up, unless the user passes in a bad
            //list of vertices. . . (Same goes for returning VERTICAL_SLOPE in
            //edgeSlope() above.)
            return VERTICAL_SLOPE;
        }
        return (double)(v2.second - v1.second) / (double)(v2.first - v1.first);
    }
#endif



//CIPs
CIP pathFromVertices(const vector<Vertex>& verts)
{
    //If we have no vertices, then we don't even have the trivial path, so throw an exception.
    if (verts.size() == 0)
        throw invalid_argument("No vertices specified.");
    //First, check that the vertices start on the y-axis and end on the x-axis.
    if (verts[0].first != 0 || verts[0].second < 0)
        throw invalid_argument("First vertex is not on y-axis.");
    if (verts.back().first < 0 || verts.back().second != 0)
        throw invalid_argument("Last vertex is not on x-axis.");
    //If there is only one vertex, then the above can only have not thrown an exception if the
    //vertex is (0,0), so that this is the trivial path. Since vector<Edge> will initialize any elements
    //to (0,0) if we use its 1-parameter construction, we can easily construct this path.
    if (verts.size() == 1) { return CIP(1); }

    /*The last things to check about the vertices are that the x-coordinates are non-decreasing and that
      the slopes of the edges between the vertices are nonpositive and non-increasing. If these are true,
      then the edges form a valid CIP, which we will go ahead and create as we check these conditions.*/
    CIP path;
    /*The maximum number of edges we could have is the number of vertices minus 1 (note that there will be
      less than this if more than two vertices in a row are on the same edge). So, we reserve space for that
      many edges to make sure we don't have to resize.*/
    path.reserve(verts.size() - 1);
    double lastslope = 0;
    double thisslope;
    Edge nextedge;
    for (size_t i = 0, topi = verts.size() - 1; i < topi; i++)
    {
        //Compute the next edge, and make sure its x-coordinate is nonnegative (so that the x-coordinates
        //of the vertices are non-decreasing).
        nextedge = verts[i + 1] - verts[i];
        if (nextedge.first < 0)
            throw invalid_argument("x-coordinates of vertices are not non-decreasing.");
        /*Make sure the slopes are still non-increasing. Because lastslope is initialized to 0, the first
          iteration will check that the first edge has has nonpositive slope. As long as the slopes are
          non-increasing, this will imply that all the slopes are nonpositive, so we can skip checking the
          sign of the slope separately here.*/
        thisslope = edgeSlope(nextedge);
        if (thisslope > lastslope)
            throw invalid_argument("Slopes of edges are not either not nonnegative or not non-increasing.");

        //If this edge has a different slope than the last one, add it to the CIP; otherwise, increase the
        //multiplicity of the last edge by adding this edge to it. Of course, on the first iteration, there
        //is no previous edge, so we don't have to worry about combining edges.
        if (i != 0 && thisslope == lastslope)
        {
            //TODO: += is not defined here, so this seems to be constructing a new Edge every time
            //for no reason. Should we define our own class for Edge so that we can define +=?
            //(The same issue occurs in pathFromEdges().)
            path.back() = path.back() + nextedge;
        }
        else
            path.push_back(nextedge);
        lastslope = thisslope;
    }
    //At this point, we know that our path is a valid one. Before we return it, though, we wil shrink the vector
    //to fit the edges in it. This usually makes sense, since we don't expect CIP's to change after they're created.
    path.shrink_to_fit();
    return path;
}
CIP pathFromEdges(const vector<Edge>& edges)
{
    //If we have no edges, then we don't even have the trivial path, so throw an exception.
    if (edges.size() == 0)
        throw invalid_argument("No edges specified.");
    //Our loop here is essentially the same as that in pathFromVertices() above, except that we
    //don't need to compute the edges. We do still need to construct a new CIP, though, because consecutive
    //edges in the input vector might have the same slope, in which case we should combine them into one edge.
    CIP path;
    path.reserve(edges.size());
    double lastslope = 0;
    double thisslope;
    for (size_t i = 0, topi = edges.size(); i < topi; i++)
    {
        thisslope = edgeSlope(edges[i]);
        if (edges[i].first < 0)
            throw invalid_argument("Edges do not have nonnegative x-coordinate.");
        if (edges[i].second > 0)
            throw invalid_argument("Edges do not have nonpositive y-coordinate.");
        if (thisslope > lastslope)
            throw invalid_argument("Slopes of edges are not either not nonnegative or not non-increasing.");
        if (i != 0 && thisslope == lastslope)
            path.back() = path.back() + edges[i];
        else
            path.push_back(edges[i]);
        lastslope = thisslope;
    }
    //At this point, we know that our path is a valid one. Before we return it, though, we wil shrink the vector
    //to fit the edges in it. This usually makes sense, since we don't expect CIP's to change after they're created.
    path.shrink_to_fit();
    return path;
}

//CGen constructors
CGen::CGen() : path(1), xval(0), yval(0), mval(0), hval(0)
{}
CGen::CGen(const CIP& p, const vector<char>& ls) : path(p), labels(ls)
{
    //Because CGen's are not typically changed (we just construct new ones), we may as
    //well make sure our path and labels are not holding extra memory.
    path.shrink_to_fit();
    labels.shrink_to_fit();
    xval = 0; yval = 0; mval = 0; hval = 0;
    for (size_t i = 0, topi = p.size(); i < topi; i++)
    {
        xval += p[i].first;
        yval -= p[i].second;
        mval += boost::integer::gcd(p[i].first, p[i].second);
        if (ls[i] != 0)
            hval++;
    }
}
CGen CGenFromVertices(const vector<Vertex>& verts, const vector<char>& labels)
{
    //First, use pathFromVertices() to check if the vertices are valid and create them if so.
    CIP path = pathFromVertices(verts);
    //Now, we need to make sure that we have the right number of labels and that horizontal and
    //vertical edges are not labelled "h."
    if (path.size() != labels.size())
        throw invalid_argument("Number of labels is not equal to the number of edges.");
    //Notice that the only the first edge can be horizontal (i.e. have slope 0), and only the last edge can be
    //vertical (i.e. have slope -infinity).
    if (path[0].second == 0 && labels[0] == 'h')
        throw invalid_argument("Horizontal edge is labelled 'h.'");
    if (path.back().first == 0 && labels.back() == 'h')
        throw invalid_argument("Vertical edge is labelled 'h.'");
    //The final thing we have to check is that all of our edges are labelled either 'e' or 'h' and not something
    //else entirely. We also need to create a vector that represents labels as 0's and 1's rather than as the
    //characters 'e' and 'h,' so we will do this while we check that the labels are correct.
    vector<char> newlabels;
    newlabels.reserve(path.size());
    for (size_t i = 0, topi = path.size(); i < topi; i++)
    {
        if (labels[i] == 'e')
            newlabels.push_back(0);
        else if (labels[i] == 'h')
            newlabels.push_back(1);
        else
            throw invalid_argument("Edge is not labelled 'e' or 'h.'");
    }
    return CGen(path, newlabels);
}
CGen CGenFromEdges(const vector<Edge>& edges, const vector<char>& labels)
{
    //First, use pathFromEdges() to check if the edges are valid and create a CIP if so.
    CIP path = pathFromEdges(edges);
    //Now, we just check that the labels are valid and reformat them. This is identical to
    //what we do in CGenFromVertices().
    if (path.size() != labels.size())
        throw invalid_argument("Number of labels is not equal to the number of edges.");
    if (path[0].second == 0 && labels[0] == 'h')
        throw invalid_argument("Horizontal edge is labelled 'h.'");
    if (path.back().first == 0 && labels.back() == 'h')
        throw invalid_argument("Vertical edge is labelled 'h.'");
    vector<char> newlabels;
    newlabels.reserve(path.size());
    for (size_t i = 0, topi = path.size(); i < topi; i++)
    {
        if (labels[i] == 'e')
            newlabels.push_back(0);
        else if (labels[i] == 'h')
            newlabels.push_back(1);
        else
            throw invalid_argument("Edge is not labelled 'e' or 'h.'");
    }
    return CGen(path, newlabels);
}
//A helper function for parseCGen() below.
void addEdge(vector<Edge>* edges, vector<char>* labels, boost::fusion::vector4<int, char, int, int> edge)
{
    using boost::fusion::at_c;
    edges->push_back(make_pair(at_c<0>(edge) * at_c<2>(edge), at_c<0>(edge) * (-abs(at_c<3>(edge))) ));
    labels->push_back(at_c<1>(edge));
}
template <typename Iterator>
bool parseCGen(const Iterator& first, const Iterator& last, CGen& target)
{
    //Our normal parsing will not work if we are looking for the trivial CGen, which is represented
    //by "1." So, we handle this case separately.
    if (string(first, last) == "1")
    {
        target = CGen(CIP(1), vector<char>());
        return true;
    }
    //First, some using and namespace directives to make our Spirit calls more readable.
    namespace qi = boost::spirit::qi;
    using boost::spirit::ascii::char_;
    using qi::int_;
    using qi::attr;
    using qi::parse;

    //We will need to store all of the edges and labels as we parse.
    vector<Edge> es;
    vector<Edge>* edges = &es;
    vector<char> ls;
    vector<char>* labels = &ls;
    //This is the rule that we will use to parse individual edges.
    qi::rule<Iterator, boost::fusion::vector4<int, char, int, int>()> edge =
        (int_ | attr(1)) >> char_ >> '(' >> int_ >> ',' >> int_ >> ')';
    //Now, the parse statement.
    bool success = parse(first, last, +(edge[boost::bind(addEdge,edges,labels,boost::placeholders::_1)]));

    //If we didn't manage to parse the input successfull, return false.
    if (!success) { return false; }
    //Once we finish parsing, we need to make sure the data we got is valid and then make a CGen object.
    //CGenFromEdges() will take care of this for us. Note that this will throw an std::invalid_argument exception
    //if the CGen is not valid.
    target = CGenFromEdges(es, ls);
    return true;
}

CGen CGenFromString(const string& str)
{
    CGen cg;
    const auto& first = str.begin();
    const auto& last = str.end();
    bool success = parseCGen(first, last, cg);
    if (!success)
        throw invalid_argument("Invalid string passed in");
    return cg;
}

//CGen member functions
CIP CGen::getPath() const { return path; }
vector<char> CGen::getLabels() const { return labels; }
int CGen::x() const { return xval; }
int CGen::y() const { return yval; }
int CGen::m() const { return mval; }
int CGen::h() const { return hval; }
int CGen::L() const
{
    /*We will use Pick's Theorem to compute the number of lattice points enclosed in the polygon
      bounded by our path. The statement of the theorem is that
      A = i+b/2-1,
      where A is the area of the polygon, i is the number of lattice points in the interior of the
      polygon, and b is the number of lattice points on the boundary of the polygon. Rearranging this
      equation and adding b/2 to both sides tells us that
      L = i+b = A+b/2+1.
      So as not to deal with fractions, we will compute 2A+b (which is necessarily even because the above
      equation must give us an integer) and divide this by 2 to get the final result.

      Now, the number of boundary points is x()+y()+m(). To compute the area (or rather, 2 times the area),
      we will sum up the area under each of the edges.*/
    int A2 = 0;
    //We need to know how high up each edge is in order to know its area. For this, we will iterate backwards through
    //the edges. The last edge must end on the x-axis, which tells us its height; the rest we can find by adding the
    //y-coordinates of the edges we've iterated through already.
    int ysum = 0;
    for (auto it = path.rbegin(); it != path.rend(); ++it)
    {
        //If this edge is vertical, there is no area underneath it.
        if (it->first != 0)
        {
            /*The part of the polygon underneath this edge is a trapezoid. (This isn't quite true -- if we are on
              the last edge, then it will be a triangle, but the way that the area formulas work out, we will get the
              same answer by regarding a triangle as a trapezoid with one of the parallel sides having length 0. It
              could also be a rectangle, but this is simply a special kind of trapezoid, so again, the same area formula
              works.) So, the formula for twice its area is
              x*(y1+y2),
              where x is the width of the trapezoid (i.e. the x-coordinate of the edge), y1 is the y-coordinate of
              the upper-left vertex of the edge, and y2 is the y-coordinate for the lower-right vertex of the edge.
              Here, y2 is simply ysum, while y1 is ysum - the y-coordinate of the edge (minus because the y-coordinate
              is negative).*/
            A2 += it->first * (2 * ysum - it->second);
        }
        //Update ysum by adding the difference in y from this edge to it.
        ysum -= it->second;
    }
    //Now that we know what twice the area is, we can compute L immediately using our earlier comments.
    return (A2 + xval + yval + mval) / 2 + 1;
}
int CGen::I() const { return 2*(L() - 1) - hval; }
string CGen::toString() const
{
    //We will build up our string by using a stringstream, which will take care
    //of integer-to-string conversions for us.
    stringstream ss;
    int mult;
    for (size_t i = 0, topi = path.size(); i < topi; i++)
    {
        mult = boost::integer::gcd(path[i].first, path[i].second);
        //We adopt the convention of not writing the multiplicity when it is 1. In this case, we also don't
        //need to divide the edge coordinates by the multiplicity, since doing so would not change them.
        if (mult == 1)
            ss << (labels[i] == 0 ? 'e' : 'h') << '(' << path[i].first << ',' << path[i].second << ')';
        else
            ss << mult << (labels[i] == 0 ? 'e' : 'h') << '(' << (path[i].first / mult) << ','
                << (path[i].second / mult) << ')';
    }
    return ss.str();
}
CGen CGen::operator*(const CGen& rhs)
{
    /*We loop over the edges of both the input CGens, adding the edge with the next lowest slope
      (since the paths are sorted, this is at the beginning of one of the two CGens) to the resultant
      path each time and combining labels until we run out of one of the two paths. After that, we
      copy the rest of the remaining path and labels into the resultant path and labels.

      We initialize the resulting path and label vectors to be the sum of the sizes of the original
      paths. This is the maximum size they could be; they will be smaller if and only if edges of the
      same slope appear in both cg1 and cg2.*/
    int len1 = this->path.size();
    int len2 = rhs.path.size();
    vector<Edge> respath;
    respath.reserve(len1 + len2);
    vector<char> reslabels;
    reslabels.reserve(len1 + len2);
    int i1 = 0;
    //Stop whenever we run out of edges of cg2.
    for (int i2 = 0; i2 < len2; )
    {
        //Also stop whenever we run out of edges of cg1.
        if (i1 == len1)
        {
            std::copy(rhs.path.begin()+i2, rhs.path.end(), back_inserter(respath));
            std::copy(rhs.labels.begin() + i2, rhs.labels.end(), back_inserter(reslabels));
            return CGen(respath, reslabels);
        }
        //If the next edge to add is shared by both paths, then it must be the case that this
        //edge is next in both of the paths. So, we can simply compare the slopes of the next
        //two edges.
        if (edgeSlope(this->path[i1]) == edgeSlope(rhs.path[i2]))
        {
            //If both edges have the same slope, we simply add them to get the combined edge
            //that goes in the resultant path.
            respath.push_back(this->path[i1] + rhs.path[i2]);
            /*If both paths have this edge labelled 'h,' then the product is not defined, so we
              throw an exception. If one of the paths has this edge labelled 'h,' then the resultant
              path will as well; and if both paths have this edge labelled 'e,' then so will the
              resultant path.*/
            if (this->labels[i1] != 0)
            {
                if (rhs.labels[i2] != 0)
                    throw invalid_argument("The input CGens share a hyperbolic orbit.");
                reslabels.push_back(1);
            }
            else if (rhs.labels[i2] != 0)
                reslabels.push_back(1);
            else
                reslabels.push_back(0);
            //Finally, advance both i1 and i2 past this edge.
            i1++; i2++;
        }
        else if (edgeSlope(this->path[i1]) > edgeSlope(rhs.path[i2]))
        {
            //In this case, we just have to add the edge and label from this path to the resultant path.
            respath.push_back(this->path[i1]);
            reslabels.push_back(this->labels[i1]);
            //Advance i1 path the edge we just added.
            i1++;
        }
        else //edgeSlope(this->path[i1] < edgeSlope(rhs.path[i2]))
        {
            //In this case, we just have to add the edge and label from the rhs to the resultant path.
            respath.push_back(rhs.path[i2]);
            reslabels.push_back(rhs.labels[i2]);
            //Advance i2 path the edge we just added.
            i2++;
        }
    }
    std::copy(this->path.begin()+i1, this->path.end(), back_inserter(respath));
    std::copy(this->labels.begin() + i1, this->labels.end(), back_inserter(reslabels));
    return CGen(respath, reslabels);
}

//CTDs
int CTD::LT_boundX(const int index, const double action, const int hcond) const
{
    /*To put a bound on x() and y() that will work for any CTD, we need to use the index condition
      so as to avoid using the action condition of less than. Given any CGen lt such that lt <= cg,
      we must have
      lt.I() = cg.I(). Now, suppose that
      lt.y() >= cg.I()/2 + 1 = lt.I()/2 + 1.
      Then, lt.L() > lt.y() (the edge on the y-axis alone contributes lt.y() + 1 lattice points), so
      we have
      lt.I() = 2*(lt.L() - 1) - lt.h()/2 > 2*(lt.y() - 1) - lt.h()/2 >= 2*(lt.I()/2 + 1 - 1) - lt.h()/2
      >= lt.I().
      But this last equation states that lt.I() > lt.I(), which is clearly impossible. This proves that
      lt.y() < cg.I()/2 + 1,
      so that the maximum lt.y() can be is (the floor of) cg.I()/2. The exact same argument goes through for
      lt.x(), so we have maximum values for both lt.x() and lt.y().

      As discussed above, the value we want to return here is the floor of index/2. The easiest way to
      compute this is just index >> 1.*/
    int indexbound = index >> 1;
    //We also get our bound using the action from LT_boundXUsingAction().
    int actionbound = LT_boundXUsingAction(index, action, hcond);
    return min(indexbound, actionbound);
}
int CTD::LT_boundY(const int index, const double action, const int hcond) const
{
    //See comments in LT_boundX() above for explanation of this bound.
    int indexbound = index >> 1;
    //We also get our bound using the action from LT_boundYUsingAction().
    int actionbound = LT_boundYUsingAction(index, action, hcond);
    return min(indexbound, actionbound);
}
int CTD::LT_boundXUsingAction(const int index, const double, const int) const
{
    /*There isn't currently a default implementation of this method. Perhaps it
      should just be abstract; but just in case someone doesn't want to
      implement their own override, we set it to return the same index bound as
      the one in LT_boundX(). That way, not overriding this method will make it
      have no effect.*/
    return index >> 1;
}
int CTD::LT_boundYUsingAction(const int index, const double, const int) const
{
    //See LT_boundXUsingAction() for comments.
    return index >> 1;
}
Polydisk::Polydisk(double x, double y) : a(x), b(y)
{
    //Although we have already assigned x and y to a and b, we have to make sure that these
    //doubles are positive.
    if (x <= 0)
        throw invalid_argument("Polydisk has nonpositive x-value");
    if (y <= 0)
        throw invalid_argument("Polydisk has nonpositive y-value");
}
double Polydisk::actionOf(const CGen& cg) const
{
    //The action with respect to a Polydisk has a straightforward formula.
    return b*cg.x() + a*cg.y();
}
string Polydisk::toString() const
{
    return "Polydisk_" + to_string(a) + "_" + to_string(b);
}
int Polydisk::LT_boundXUsingAction(const int, const double action, const int hcond) const
{
    /*We have a relatively simple expression for actions of CGens with respect to Polydisks, and we
      can use this to get pretty good upper bounds on x() and y(). Suppose ctd1 is the polydisk P(a,b),
      and let lt be a CGen such that lt <= cg. Then, the action condition of lt <= cg tells us that
      (1) b*lt.x() + a*lt.y() <= action,
      where action is the action of cg with respect to ctd2 (i.e. it is the parameter action passed to
      both xbound() and ybound()). On the other hand, the h condition of lt <= cg tells us that
      (2) lt.x() + lt.y() >= lt.x() + lt.y() - lt.h()/2 >= hcond,
      where hcond is cg.x() + cg.y() + cg.m() - 1 (i.e. it is the parameter hcond passed to
      xbound() and ybound()).

      Suppose first that a > b. Then, (2) tells us that
      lt.x() >= hcond - lt.y().
      Substituting this into (1) gives
      b*(hcond - lt.y()) + a*lt.y() <= b*lt.x() + a*lt.y() <= action,
      or, rearranging:
      (a - b)*lt.y() <= action - b*hcond.
      Dividing by the (nonnegative) value a - b then gives us an upper bound on lt.y(). We can then get a
      (slightly weaker) bound for lt.x() straight from (1):
      b*lt.x() <= b*lt.x() + a*lt.y() <= action,
      so we can divide by b to get an upper bound on lt.x().

      Now, the same argument goes through when a < b if we swap a with b and lt.x() with lt.y(). The only
      remaining case is when a = b. In this case, our above argument for the bound on lt.y() does not go
      through, since the bound we get becomes
      0 <= action - b*hcond,
      which does not yield any bound on lt.y(). Instead, we must use our above argument for a bound on
      lt.x() to get bounds for both x() and y(). This yields
      b*lt.x() <= action
      and
      a*lt.x() <= action.
      
      One final detail: mathematically speaking, we should be able to apply floor to our upper bound 
      and maintain a valid bound (because x() and y() are integers). However, floating-point arithmetic
      errors can sometimes make floor() return the wrong value. It's better to return a weaker bound
      than an inaccurate one, so we use round() instead to avoid this problem (with the assumption
      that even all the accumulated floating-point error in the action computation is unlikely to make
      round() give an incorrect bound).
      
      For an explicit example in which this floating-point error mattered: when
      finding less thans for 2e(1,-1) with respect to Polydisk(2.43, 1) and Ellipsoid(3.215, 3.215),
      my tests failed to find the less than e(4,-1) because Polydisk::LT_boundYUsingAction() was
      rounding down its return value to 0 (the actual value should have been exactly 1, but
      with floating-point error, this became slightly less than 1, so floor made it 0).
      Changing from floor() to round() solved the problem in this case.
    */
    if (b > a) { return round((action - a*hcond) / (b - a)); }
    else { return round(action / b); }
}
int Polydisk::LT_boundYUsingAction(const int, const double action, const int hcond) const
{
    //See the comments in LT_boundXUsingAction() for explanation of this bound.
    if (a > b) { return round((action - b*hcond) / (a - b)); }
    else { return round(action / a); }
}
double Polydisk::LT_actionOf(const CIP&, const int_fast16_t, const int_fast16_t,
        const int_fast16_t x, const int_fast16_t y) const
{
    //See comments in Polydisk::actionOf().
    return b*x + a*y;
}
Ellipsoid::Ellipsoid(double x, double y) : a(x), b(y)
{
    //Although we have already assigned x and y to a and b, we have to make sure that these
    //doubles are positive.
    if (x <= 0)
        throw invalid_argument("Ellipsoid has nonpositive x-value");
    if (y <= 0)
        throw invalid_argument("Ellipsoid has nonpositive y-value");
    //Finally, we store the slope of the line that defines the Ellipsoid as a CTD.
    slope = -b / a;
}
double Ellipsoid::actionOf(const CGen& cg) const
{
    /*The formula for the action with respect to an ellipsoid is
      b*x+a*y,
      where (x,y) is a point on the line tangent to the CGen of slope -b/a. Such a tangent line goes
      through the upper left vertex of the first edge of the CGen with slope less than -b/a, or through
      the last vertex of the CGen if no such edge exists; we will use this vertex as (x,y). To find it, we
      iterate through each edge of the CGen, adding up x- and y-coordinates as we go in order to calculate
      the coordinates of the desired vertex; we stop whenever we hit an edge of slope at least as steep as
      (i.e. less than or equal to) -b/a.*/
    int xsum = 0;
    //To find the y-coordinate of the vertex, we start with the total y() of the CGen and subtract off
    //y-coordinates of edges as we go.
    int ysum = cg.y();
    //First, check if no edge has slope less than -b/a. If this is the case, we are looking for the vertex
    //(cg.x(),0), so we don't have to iterate over the edges.
    if (slope < edgeSlope(cg.getPath().back()))
        return b*cg.x();
    //Now, iterate over the edges (other than the last one, since we know from the above if statement that
    //we'll find a vertex before the last one).
    for (size_t i = 0, topi = cg.getPath().size() - 1; i < topi; i++)
    {
        if (slope >= edgeSlope(cg.getPath()[i]))
            return b*xsum + a*ysum;
        xsum += cg.getPath()[i].first;
        ysum += cg.getPath()[i].second;
    }
    //If we haven't returned yet, then slope is steeper than all the edges of the CGen except the last one.
    //So, the vertex that we want is the second-to-last one, which is currently (xsum,ysum).
    return b*xsum + a*ysum;
}
string Ellipsoid::toString() const
{
    return "Ellipsoid_" + to_string(a) + "_" + to_string(b);
}
int Ellipsoid::LT_boundXUsingAction(const int, const double action, const int) const
{
    /*We have a succinct but not entirely simple expression for the action of a CGen with respect to
      an Ellipsoid, and we can use this to find decent bounds on the x() and y() of CGens in the less than
      searcg. To see this, suppose that ctd1 is the ellipsoid E(a,b), and let lt be a CGen such that
      lt <= cg. Then, the action condition of lt <= cg tells us that
      (1) b*x + a*y <= action,
      where action is the parameter passed to xbound() and ybound() and (x,y) is the vertex of lt at which a
      line of slope -b/a is tangent. Now, for any vertex (m,n) of lt, we claim that b*m + a*n <= b*x + a*y,
      with equality if and only if (m,n) = (x,y). To see this, note that the tangent line of slope -b/a at
      (x,y) has equation

      b*w + a*z = c,

      where w and z are variables and c is some constant. Since this line is tangent to the CGen lt, we
      know that every vertex of lt must lie on one side of this line; since CGen's are convex, this must
      be the "lower left side," i.e. the half-plane defined by

      b*w + a*z <= c.

      In particular, plugging in (m,n) gives us

      b*m + a*n <= c = b*x + a*y.

      (The equality on the right here follows from the fact that (x,y) is on the tangent line by definition.)
      This inequality is an equality if and only if (m,n) is on the tangent line; but the only point in
      lt that is on the tangent line is (x,y), so equality occurs if and only if if (m,n) = (x,y).
      
      Now, applying the above when (m,n) is the vertex (lt.x(),0) and using (1) gives
      b*lt.x() + 0 <= b*x + a*y <= action,
      which gives us an upper bound on lt.x(). Likewise, taking (m,n) to be the vertex (0,lt.y()) yields
      0 + a*lt.y() <= action,
      so we also get an upper bound on lt.y().
      
      One final detail: mathematically speaking, we should be able to apply floor to our upper bound 
      an maintain a valid bound (because x() and y() are integers). However, floating-point arithmetic
      errors can sometimes make floor() return the wrong value. It's better to return a weaker bound
      than an inaccurate one, so we use round() instead to avoid this problem (with the assumption
      that even all the accumulated floating-point error in the action computation is unlikely to make
      round() give an incorrect bound). 
      */
    return round(action / b); 
}
int Ellipsoid::LT_boundYUsingAction(const int, const double action, const int) const
{
    //See the comments in Ellipsoid::LT_boundXUsingAction() for an explanation of this bound.
    return round(action / a);
}
double Ellipsoid::LT_actionOf(const CIP& path, const int_fast16_t lastind, const int_fast16_t,
        const int_fast16_t x, const int_fast16_t y) const
{
    /*See comments in Ellipsoid::actionOf(). Note that even if a final vertical
      edge is missing from path, it doesn't change anything about these
      computations: the slope -b/a is always greater than the slope of a
      vertical edge, so throwing out the vertical edge doesn't affect where the
      line of slope -b/a is tangent to the path.*/
    int xsum = 0, ysum = y;
    if (slope < edgeSlope(path[lastind]))
        return b*x;
    for (int_fast16_t i = 0; i < lastind; i++)
    {
        if (slope >= edgeSlope(path[i]))
            return b*xsum + a*ysum;
        xsum += path[i].first;
        ysum += path[i].second;
    }
    return b*xsum + a*ysum;
    return 0;
}
LinearCTD::LinearCTD(const CIP& path)
{
    /*We mainly need to compute the vertices that lie on the given path. For this, we can sum up
      the x- and y-coordinates of the edges. The only issue is that we need to sum edges from
      left-to-right to get x-coordinates and from right-to-left to get y-coordinates. We can either
      loop once and create two different arrays of coordinates or loop twice and create one array of
      coordinates. Since performance shouldn't really matter in this function, we'll do the latter, because
      it's more straightforward.*/
    vector<int> ys;
    //We will necessarily have one more vertex than we have edges. However, we don't need to save the
    //last y-value, since it's always 0, so we need exactly as many y-values as we have edges.
    ys.reserve(path.size());
    int total = 0;
    for (auto it = path.rbegin(); it != path.rend(); ++it)
    {
        total -= it->second;
        ys.push_back(total);
    }

    //Now that we have the y-values, we can compute the x-values and create the vertices at once. We will
    //also go ahead and store the slopes of all the edges.
    vertices.reserve(path.size() + 1);
    slopes.reserve(path.size());
    total = 0;
    //We don't do the lasat edge yet, because we will not save its slope or
    //final vertex if it is vertical.
    for (size_t i = 0, topi = path.size() - 1; i < topi; i++)
    {
        //Because we created the y-values in reverse order, we need to access them in reverse order here.
        vertices.push_back(make_pair(total, ys[ys.size() - 1 - i]));
        slopes.push_back(edgeSlope(path[i]));
        total += path[i].first;
    }
    //Store the second-to-last vertex (the upper left vertex of the last edge).
    vertices.push_back(make_pair(total, ys[0]));
    //Now, we've stored everything except the last vertex and the last slope. We
    //store these if and only if the last edge is not vertical.
    if (path.back().first != 0)
    {
        vertices.push_back(make_pair(total + path.back().first, 0));
        slopes.push_back(edgeSlope(path.back()));
    }
}
LinearCTD CTDFromVertices(const vector<Vertex>& verts) { return LinearCTD(pathFromVertices(verts)); }
LinearCTD CTDFromEdges(const vector<Edge>& edges) { return LinearCTD(pathFromEdges(edges)); }
double LinearCTD::actionOf(const CGen& cg) const
{
    /*The general formula for the action of a CTD on a CGen is the sum over every edge (x,y) of the
      CGen of x*b-y*a, where (a,b) is a point on the line tangent to the CTD and parallel to the
      edge (x,y). So, we simply iterate over each edge of the CGen to sum all such terms.*/
    double action = 0;
    size_t lastvert = 0;
    CIP path = cg.getPath();
    double eslope;
    for (auto it = path.begin(); it != path.end(); ++it)
    {
        eslope = edgeSlope(*it);
        //Because the edges of the CGen are in order of decreasing slope, each one will be tangent at a
        //vertex to the right of (or equal to) the vertex at which the previous edge was tangent. This
        //allows us to initialize our vertex index at the beginning and not reset it for different edges.
        for (; lastvert < vertices.size() - 1 && slopes[lastvert] >= eslope; lastvert++) { }
        /*Now, lastvert is either at the vertex we want or at the last vertex. In the latter case, the CGen
          edge has a slope at least as big as all the edges of the CTD boundary (except possibly a final
          vertical one), so the line with the edge slope is tangent at the last vertex. Either way,
          lastvert is now the vertex we want.*/
        action += (*it).first * vertices[lastvert].second - (*it).second * vertices[lastvert].first;
    }
    return action;
}
string LinearCTD::toString() const
{
    stringstream ss("LinearCTD_");
    for (const Vertex& v : vertices)
        ss << "(" << v.first << "," << v.second << ")-";
    return ss.str();
}
double LinearCTD::LT_actionOf(const CIP& path, const int_fast16_t lastind, const int_fast16_t vertedge,
        const int_fast16_t, const int_fast16_t) const
{
    //See comments in LinearCTD::actionOf().
    double action = 0;
    size_t lastvert = 0;
    double eslope;
    for (int_fast16_t i = 0; i <= lastind; i++)
    {
        eslope = edgeSlope(path[i]);
        for (; lastvert < vertices.size() - 1 && slopes[lastvert] >= eslope; lastvert++) {}
        action += path[i].first * vertices[lastvert].second - path[i].second * vertices[lastvert].first;
    }
    /* If a vertical edge is present in the CGen, then it isn't in path, so we haven't dealt with it yet.
       Such an edge is always tangent at the final vertex of the path that bounds the CTD. In this case,
       the vertical edge is given by (0, -vertedge), so the contribution of this vertical edge to the action is:
       0 * vertices.back().second - (-vertedge) * vertices.back().first = vertedge * vertices.back().first.*/
    if (vertedge != 0)
        action += vertedge * vertices.back().first;
    return action;
}
