//data.h : Here we define the main data types we will use. The bulk of the API of the program
//lies in two of these data types, namely CGen and CTD. However, there are also some more fundamental
//data types in this file.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

using namespace std;

/************************************************************************************************************
  Vertices and Edges
 ***********************************************************************************************************
 The most fundamental data types we will work with are vertices and edges. Vertices are represented
 by their x- and y-coordinates, while edges are represented by the x- and the y-coordinates of the
 vector pointing from the start of the edge to the end of the edge. Note from the type definition that
 we represent our coordinates with the type int_fast16_t. This means that we make no guarantees
 on the exact size of these numbers except that they are signed and at least 16 bits. To avoid overflow,
 then, make sure you don't work with any coordinates larger than 32,767.*/
typedef pair<int_fast16_t, int_fast16_t> Vertex, Edge;
//Here we define component-wise addition on edges, so that we can treat them like 2-dimensional vectors.
Edge operator+(const Edge&, const Edge&);
Edge operator-(const Edge&, const Edge&);
//Compute the slope of an edge, i.e. its y-coordinate divided by its x-coordinate.
inline double edgeSlope(const Edge&);
//Compute the slope of the line segment between two vertices.
inline double vertexSlope(const Vertex&, const Vertex&);

/************************************************************************************************************
  Convex Integral Paths (CIPs)
 ***********************************************************************************************************
 A convex integral path (CIP) is the boundary of a convex polygon P whose vertices have integer coordinates
 such that
 1. for some nonnegative integers x and y, the line segments from (0,0) to (x,0) and from (0,0) to (0,y)
 are edges of P, and
 2. P lies inside the rectangle with vertices (0,0), (x,0), (0,y), and (x,y).
 (This is not the original definition but rather a more geometric one that I believe to be equivalent.
 For the original, see "Beyond ECH Capacities," p. 5, Definition 1.9.)
 We will represent such a polygon by a list of all of its edges except for the one from (0,0) to (x,0) and
 the one from (0,0) to (0,y). Notice that we can recover x (respectively y) from such a list by taking the
 sum of all the x-components (respectively y-components) of the edges in the list.

 Now, it turns out that we are mainly interested in a specific kind of CIP, namely a convex generator. We define
 a class to represent this object below. Because we don't plan to work with CIP's alone very much, we will simply
 give them a type definition as a vector of edges and not enforce that the edges in question define a
 valid CIP. We do, however, provide two functions, pathFromEdges() and pathFromVertices(), that create CIP's and
 enforce proper type safety for them. If you are looking to use CIP's alone, you should use those functions.*/
typedef vector<Edge> CIP;
//These two functions create CIP's from vectors of vertices and edges, respectively. These return a valid CIP if
//their inputs are valid and throw an std::invalid_argument exception otherwise.
CIP pathFromVertices(const vector<Vertex>&);
CIP pathFromEdges(const vector<Edge>&);

/************************************************************************************************************
  Convex Generators (CGens)
 ***********************************************************************************************************
 This class represents a convex generator. A convex generator is a CIP along with a label of either "e" or
 "h" for each edge such that no horizontal or vertical edges are labelled "h."*/
class CGen
{
    private:
        /*In addition to storing the CIP corresponding to the CGen, we store a vector where each element
          is 0 if the edge at the same position in the CIP is labelled "e" and a 1 if the edge is labelled "h."
          Note that we are using vector <char> here rather than vector<bool>: the latter is controversial, so
          this is a standard workaround when space is not a crucial consideration (as in this case, since we don't
          expect to have a very large number of CGens or a very large number of edges per CGen).
          We also store the values of x(), y(), m(), and h(): we use these to compute L() and I(), so they are
          frequently useful, and since computing them requires iterating through the whole CIP, we would rather
          not recalculate them every time we need them.*/
        CIP path;
        vector<char> labels;
        int xval, yval, mval, hval;
        /*A private constructor that simply sets path and labels to the specified arguments and then calculates
          x, y, m, and h. This is used by CGenFromEdges, CGenFromVertices, and CGenFromString. Note that this
          constructor does NOT check that its inputs define a valid CGen; that functionality is handled by the
          above-mentioned methods.*/
        CGen(const CIP&, const vector<char>&);
    public:
        //This default constructor creates the trivial CGen.
        CGen();
        /*These friend functions create CGens from given vertices and edges, respectively. They use
          pathFromVertices() and pathFromEedges(), respectively, to check that their first parameters are valid
          and to create a CIP if so. The second parameter must contain an 'e' at each index corresponding to an
          edge labelled 'e' and an 'h' at each index corresponding to an edge labelled 'h'; no other char values
          can be used. Naturally, this parameter must have the same length as the number of edges of the desired
          path (which is 1 less than the number of vertices of the path, i.e. points at which the slope changes).
          Moreover, we cannot have a label of "h" corresponding to any horizontal or vertical edge. If both
          parameters are valid, these functions return CGen's created from the given inputs. Otherwise, they
          throw an std::invalid_argument exception.*/
        friend CGen CGenFromVertices(const vector<Vertex>&, const vector<char>&);
        friend CGen CGenFromEdges(const vector<Edge>&, const vector<char>&);
        /*This function parses a CGen from an iterator over a string representation of it. For the correct
          format of such a string, see the comment on CGen.toString() below. If the input is of the correct format,
          the last parameter is set to reference the corresponding CGen, and the function returns true; if the input
          is of the correct format but specifies an invalid CGen, then an std::invalid_argument exception is thrown,
          and the last parameter is not changed; and if the input is of an incorrect format, then the function returns
          false, and the last parameter is not changed. As is the policy of the Spirit library, the Iterator first
          is moved to the end of the section of input that can be successfully parsed. Note that if the entire input
          is not consumed, the function will still return true and change the last parameter; if this is not the desired
          behavior, the caller can compare first and last after calling the function to see if all the input was consumed.
          This function is mainly useful for parsing Cgens from various sources programmatically; for a convenient wrapper
          that parses a string into a CGen, see CGenFromString().*/
        template <typename Iterator>
            friend bool parseCGen (const Iterator&, const Iterator&, CGen&);
        //Get the CIP associated to this CGen.
        CIP getPath() const;
        //Get the labels associated to this CGen.
        vector<char> getLabels() const;
        //Get the x-value of the CGen, i.e. the x in the definition of CIP given above.
        int x() const;
        //Get the y-value of the CGen, i.e. the y in the definition of CIP given above.
        int y() const;
        //Get the total multiplicity of the edges of the CGen, which is equivalent to the
        //number of integer points on the CIP other than the origin.
        int m() const;
        //Get the total number of edges of the CGen labelled h.
        int h() const;
        //Compute the number of points with integer coordinates that lie in the polygon bounded by
        //the CIP.
        int L() const;
        /*Compute the ECH index of the CGen. this is defined as
          2(L()-1)-h(),
          where L() and h() are defined as above.*/
        int I() const;
        /*This method converts the CGen to a string that represents it. For any edge (x,y) of the CGen,
          suppose that m = gcd(x,y) is the multiplicity of the edge, a = x/m, and b = -y/m (note the use of
          the minus sign to define b, which makes both a and b are positive). Then, we represent (x,y) by
          the string "me(a,b)" if (x,y) is labelled "e" and "mh(a,b)" if (x,y) is labelled "h."
          With this convention, we can represent a CGen by making a string to represent each of its edges
          and concatenating them. We will order the edges from left to right (closest to the y-axis to farthest
          from it).*/
        string toString() const;
        /*Given any two CGens which do not share a hyperbolic orbit (i.e. which do not have edges of the same slope
          that are both labelled 'h'), then we define their product to be the CGen whose path is the path that
          contains precisely the edges of both original paths (counted with multiplicity), with each edge labelled 'h'
          if and only if that edge was labelled 'h' in one of the original CGens.*/
        CGen operator*(const CGen&);
};
//This function allows one to easily construct a CGen from a string representation of it. If the whole string
//is successfully parsed, then the corresponding CGen is returned; otherwise, an std::invalid_argument exception is thrown.
CGen CGenFromString(const string&);

/************************************************************************************************************
  Convex Toric Domains (CTDs)
 ***********************************************************************************************************
 A convex toric domain (CTD) is a region of the first quadrant bounded by the x-axis, the y-axis, a vertical
 line x = A and some non-negative, non-increasing, concave (down) function f. Since allowing f to be a general
 such function in this definition would require a lot of machinery (to take actions with respect to CTD's, we
 need to compute tangent lines, so we would have to be able to take the derivative of a general function), and
 since most of the interesting cases are very specific and nice, we make an abstract CTD class and implement
 it for a few potentially interesting cases.

 The main use of a CTD is its actionOf() method, which computes the action of a parameter CGen with respect
 to the CTD. However, we also provide a pure virtual toString() method. This method is meant to be used to
 create directories whose names are string representations of CTDs, so if you are implementing your own CTD,
 there are two main requirements that the output of toString() must satisfy. First, the returned strings
 must contain valid characters for directory names in the filesystem. Second, the output should capture all
 of the information of the CTD -- in other words, for any CTD's c1 and c2, we require that
 c1.toString() == c2.toString() if and only if c1 and c2 represent the same CTD. FSMemoizer uses toString()
 to memoize less than search results, so if you don't satisfy these requirements, then FSMemoizer may not
 work correctly for your CTD class.

 If you're just trying to use the API, that's all of the functionality you will likely need for CTD's!
 However, if you want to understand/customize how the less than search is done, this next bit will
 interest you.

 There is another, less obvious way in which we can utilize the polymorphism of CTD's. When we are trying
 to compute CGens less than a given CGen, there are 3 conditions our less than must satisfy: one condition
 on its index, one on its action, and one on its x() and y() and h() values. We would like to use these
 conditions to compute upper bounds on the possible x() and y() values our less than generator must have.
 The index condition will generally not be useful. (Pick's theorem gives us the index as a function of the
 area, x(), y(), and m(); however, because area and m() can vary widely even for fixed x() and y(), we can
 intuitively have a wide range of possible x() and y() values, all of which can be "balanced out" by a
 particular area and/or m() value to make the index condition hold.) The second and third conditions, by
 contrast, can often be used together to derive upper bounds on both x() and y(). However, to know what
 the second condition is, we must have some expression for the action of a CTD on the less than generator.
 Because this depends on the type of CTD, so do the bounds on x() and y().

 As a result, we endow the CTD base class with LT_boundX() and LT_boundY() functions that are meant to
 return the maximum possible values of x() and y() (respectively) of any CGen that satisfies the less than
 conditions specified by the function parameters. By default, these functions use a general but relatively
 weak bound on the index condition so as not to have to use any information about the action condition
 (which would depend on the specific CTD in question). However, CTD child classes can override
 LT_boundXUsingAction() and LT_boundYUsingAction() to provide additional bounds. In particular,
 Polydisk and Ellipsoid both come with their own overrides of these functions, which typically yield much
 better bounds in interesting cases.*/
class CTD
{
    protected:
        //Because CTD has child classes, it is best practice to make a protected, virtual destructor. If we
        //don't do this, then deleting a pointer to a CTD will result in undefined behavior.
        virtual ~CTD() {}
    public:
        virtual double actionOf(const CGen&) const = 0;
        virtual string toString() const = 0;
        int LT_boundX(const int index, const double action, const int hcond) const;
        int LT_boundY(const int index, const double action, const int hcond) const;
        virtual int LT_boundXUsingAction(const int index, const double action, const int hcond) const;
        virtual int LT_boundYUsingAction(const int index, const double action, const int hcond) const;
        /*This function is meant to behave just like actionOf(), but its parameters cater to the representions
          of CGens internal to LTGen in order to facilitate the less than search. The first parameter is the
          path of the CGen, possibly with extra edges on the end of it. The second parameter is the last index
          in the path that constitutes an edge that's actually part of the path; everything after this index
          should be ignored. The third parameter is the length of a vertical edge if one exists, or else 0.
          Note that if a vertical edge is part of the current result, it will NOT be on the CIP passed in, so
          one must use this parameter in order to deal with that edge properly. Finally, the fourth and fifth
          paremeters are the x() and y() of the CGen, respectively.*/
        virtual double LT_actionOf(const CIP& path, const int_fast16_t lastind, const int_fast16_t vertedge,
                const int_fast16_t x, const int_fast16_t y) const = 0;
};
/*This CTD represents a polydisk, which is a 4-dimensional generalization of a cylinder. As a CTD,
  the polydisk P(a,b) corresponds to a rectangle with vertices at (0,0), (a,0), (0,b), and (a,b).*/
class Polydisk : public CTD
{
    private:
        double a;
        double b;
    public:
        Polydisk(double, double);
        double actionOf(const CGen&) const;
        string toString() const;
        int LT_boundXUsingAction(const int index, const double action, const int hcond) const;
        int LT_boundYUsingAction(const int index, const double action, const int hcond) const;
        double LT_actionOf(const CIP& path, const int_fast16_t lastind, const int_fast16_t vertedge,
                const int_fast16_t x, const int_fast16_t y) const;
};
/*This CTD represents a (4-dimensional) Ellipsoid. As a CTD, the ellipsoid E(a,b) is a triangle with
  vertices (0,0), (a,0), and (0,b).*/
class Ellipsoid : public CTD
{
    private:
        double a;
        double b;
        //We will need the slope -b/a when computing the action with respect to the ellipsoid. Since double
        //division gets expensive over many many computations, we compute it in the constructor and
        //store it here.
        double slope;
    public:
        Ellipsoid(double, double);
        double actionOf(const CGen&) const;
        string toString() const;
        int LT_boundXUsingAction(const int index, const double action, const int hcond) const;
        int LT_boundYUsingAction(const int index, const double action, const int hcond) const;
        double LT_actionOf(const CIP& path, const int_fast16_t lastind, const int_fast16_t vertedge,
                const int_fast16_t x, const int_fast16_t y) const;
};
/*This CTD represents the case in the above definition of a CTD where f is piecewise linear, with
  all lines having rational slope and passing through at least one point with integer coordinates.
  Comparing definitions, one can see that these are precisely the CTD's that are bounded by CIP's.
  
  A few things to note about these types of CTD's (which we call LinearCTD's):
  
  (1) The corresponding CIP will have a vertical edge if and only if the function f satisfies
      f(A) > 0, where x = A is the vertical line that bounds the "rightmost edge" of the CTD.
  (2) Ellipsoids and Polydisks are actually special cases of LinearCTD's. Perhaps those
      classes should inherit from this one on principle; in practice, however, Ellipsoids
      and Polydisks have no real use for this inheritance.
  (3) The name "LinearCTD" is perhaps slightly misleading, as one could have a CTD bounded
      by a piecewise linear function f which is not bounded by a CIP (hence is not represented
      by this class). This happens if the lines that define f don't pass through points with
      integer coordinates, e.g. if f is the single line f(x) = -x + pi. */
class LinearCTD : public CTD
{
    private:
        /*Because computing actions requires vertices and not edges, we represent this CTD internally by a
          vector of vertices rather than with a CIP (which is a vector of edges). If the CIP has a vertical
          edge, we will exclude the final vertex (the one on the x-axis) from this list of vertices, because
          that vertex will never be needed for any action computations. (Actions are computed using points
          on tangent lines. Any non-vertical tangent will never be tangent at the final vertex of a
          vertical edge; vertical tangents, on the other hand, will overlap with vertical edges, so any
          point on the vertical edge can be chosen, including the second-to-last vertex. In both cases,
          the final vertex of the vertical edge is unnecessary.) */
        vector<Vertex> vertices;
        //Because edge slopes play an important role in computing actions and double division is expensive
        //if done frequently enough, we also store a vector containing the slopes of all the edges of the CTD.
        vector<double> slopes;
        /*This constructor creates a LinearCTD from a given CIP. It is used by CTDFromVertices() and
          CTDFromEdges(). Note that the constructor does NOT check if its input is a valid CIP; this is left
          to the aforementioned functions. */
        LinearCTD(const CIP& path);
    public:
        /*These two friend functions create LinearCTD's from edges and vertices, respectively. They use
          pathFromVertices() and pathFromEdges(), respectively, to create a CIP, and then they use the private
          constructor to make a LinearCTD from this path. The functions return the resulting LinearCTD if their
          input is valid and throw an std::invalid_argument exception otherwise.*/
        friend LinearCTD CTDFromVertices(const vector<Vertex>& vertices);
        friend LinearCTD CTDFromEdges(const vector<Edge>& edges);
        double actionOf(const CGen&) const;
        string toString() const;
        double LT_actionOf(const CIP& path, const int_fast16_t lastind, const int_fast16_t vertedge,
                const int_fast16_t x, const int_fast16_t y) const;
};
