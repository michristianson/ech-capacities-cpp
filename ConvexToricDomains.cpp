//ConvexToricDomains.cpp : Defines the entry point for the console application. 
#include <iostream>
#include <string>
#include <sstream>

#include "lessthans.hpp"

using namespace std;

//A little helper function to compute binomial coefficients
//(a.k.a. "n choose k"). There are likely better ways to do this, but since
//we're only using this to add up the total number of less thans we get in
//printAllLTs(), this seems good enough to me.
int choose(int n, int k)
{
    if (k == 0) { return 1; }
    return (n * choose(n-1, k-1)) / k;
}

//Functions to obtain input from stdin.
double askForDouble(string prompt)
{
    bool got_valid_input = false;
    double number;
    while (!got_valid_input)
    {
        string input;
        cout << prompt << "\n";
        getline(cin, input);
        try
        {
            number = std::stod(input);
            got_valid_input = true;
        }
        catch (invalid_argument& e)
        {
            cout << "Error: input was not a valid double." << endl;
        }
        catch (exception& e)
        {
            cout << "Error reading input: " << e.what() << endl;
        }
    }
    return number;
}
CGen askForCGen()
{
    bool got_valid_input = false;
    CGen cg;
    while (!got_valid_input)
    {
        string input;
        cout << "Please enter a string representation of a convex generator." << endl;
        getline(cin, input);
        try
        {
            cg = CGenFromString(input);
            got_valid_input = true;
        }
        catch (exception& e)
        {
            cout << "Error reading input: " << e.what() << endl;
            cout << "Some examples of valid convex generators: e(1,1)2e(1,-1)3e(0,-1) and e(4,2)h(1,3)." << endl;
        }
    }
    return cg;
}



// Functions to create nice output.
string formatCGenStats(const CGen& cg)
{
    size_t i;
    stringstream output;
    output << "CGen: " << cg.toString() << " (path = {";
    for (i = 0; i < cg.getPath().size()-1; i++)
        output << '(' << cg.getPath()[i].first << ',' << cg.getPath()[i].second << "), ";
    output << '(' << cg.getPath()[cg.getPath().size() - 1].first << ',' << cg.getPath()[cg.getPath().size() - 1].second << ')';

    output << "}, labels = {";
    for (i = 0; i < cg.getLabels().size() - 1; i++)
        output << (cg.getLabels()[i] == 0 ? 'e' : 'h') << ", ";
    output << (cg.getLabels()[cg.getLabels().size() - 1] == 0 ? 'e' : 'h');

    output << "})\n" << "statistics: x = " << cg.x() << ", y = " << cg.y() << ", m = " << cg.m()
                     << ", h = " << cg.h() << ", L = " << cg.L() << ", I = " << cg.I() << "\n";
    return output.str();
}
template<class Memoizer>
void printAllLTs(LTGen<Memoizer>& lts)
{
    int j, len;
    int i;
    int totalnum = 0;
    for (i = 1; !lts.atEnd(); i++, lts.next())
    {
        cout << i << ": Path = {";
        len = lts.currLength();
        for (j = 0; j < len - 1; j++)
            cout << '(' << lts.currPath()[j].first << ',' << lts.currPath()[j].second << "), ";
        cout << '(' << lts.currPath()[len - 1].first << ',' << lts.currPath()[len - 1].second << ')';
        cout << "}, number of h's = " << lts.currNumhs() << "\n";

        //Note: LTGen does not pick edges to label with h's. To get the total
        //number of less thans, we have to account for all possible labellings.
        int possiblehs = len;
        if (lts.currPath()[0].second == 0)
            possiblehs--;
        if (lts.currPath()[len - 1].first == 0)
            possiblehs--;
        totalnum += choose(possiblehs, lts.currNumhs());
    }
    lts.end();
    if (i == 1)
    {
        cout << "No less thans found." << endl;
    }
    else
    {
        cout << "Total number of less thans: " << totalnum << endl;
        cout << "Note: this total includes all different labelings of edges with h's." << endl;
    }
}




/*Outputs all less than results for many different computations to test our
  less than generation. There will be quite a lot of output (over 39,000 lines of
  text, which came out to a 2.4MB text file on my system), so it's probably
  best to pipe output into a file if you want to run this function. */
void testLTs()
{
    const CTD& ctd1 = Polydisk(2.43, 1);
    vector<Ellipsoid> balls{ Ellipsoid(1,1),
                       Ellipsoid(3,3),
                       Ellipsoid(3.205,3.205),
                       Ellipsoid(3.215,3.215),
                       Ellipsoid(3.225,3.225),
                       Ellipsoid(3.5,3.5),
                       Ellipsoid(4,4),
                       Ellipsoid(10,10) };
    vector<CGen> cgs;
    for (int i = 1; i <= 15; i++)
        cgs.push_back(CGenFromString(to_string(i) + "e(1,-1)"));

    for (CGen cg : cgs)
    {
        for (CTD& ctd2 : balls)
        {
            cout << "----------------------------------------------------\n"
                 << "Less thans for: " << cg.toString() << "\n"
                 << "With respect to: " << ctd1.toString() << " and " << ctd2.toString() << "\n"
                 << "----------------------------------------------------" << endl;
            LTGen<NoMemoizer> lts;
            lts.start(&ctd1, &ctd2, cg);
            printAllLTs(lts);
        }
    }

}

// Automates "less-than" computations and prints results to stdout.
void computeGivenLTs()
{
    double a = askForDouble("Please enter a value of a to define the polydisk P(a,1).");
    double c = askForDouble("Please enter a value of c to define the ball E(c,c).");
    const CTD& ctd1 = Polydisk(a, 1);
    const CTD& ctd2 = Ellipsoid(c, c);
    const CGen& cg = askForCGen();

    cout << "------------------------------------------------------\n"
         << "Info on input convex generator\n"
         << formatCGenStats(cg)
         << "------------------------------------------------------\n"
         << "Less Thans\n"
         << "----------" << endl;
    
    LTGen<NoMemoizer> lts;
    lts.start(&ctd1, &ctd2, cg);
    printAllLTs(lts);
}

int main()
{
    //testLTs();
    computeGivenLTs();
    return 0;
}
