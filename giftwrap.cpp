
// A C++ program to find convex hull of a set of points
#include <bits/stdc++.h>
#include <fstream>
using namespace std; 
ofstream fout;
struct Point					// To store the co-ordinates of every point
{
    int x, y;
} ;
// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(Point p, Point q, Point r)
{
    int val = (q.y - p.y) * (r.x - q.x) -
              (q.x - p.x) * (r.y - q.y);
    if (val == 0) return 0; 	 // colinear
    return (val > 0)? 1: 2; 	// clock or counterclock wise
}
// Prints convex hull of a set of n points
void convexHull(vector<Point> points, int n)

	{
    // There must be at least 3 points
    if (n < 3) return;
    // Initialize Result
    vector<Point> hull;
    // Find the leftmost point
    int l = 0;
    for (int i = 1; i < n; i++){
        if (points[i].x < points[l].x){
            l = i;
        }
    }
    // Start from leftmost point, keep moving counterclockwise
    // until reach the start point again.  This loop runs O(h)
    // times where h is number of points in result or output.
    int p = l, q;
    do
    {
        // Add current point to result
        hull.push_back(points[p]);
        // Search for a point 'q' such that orientation(p, x,
        // q) is counterclockwise for all points 'x'. The idea
        // is to keep track of last visited most counterclock-
        // wise point in q. If any point 'i' is more counterclock-
        // wise than q, then update q.
        q = (p+1)%n;
        for (int i = 0; i < n; i++)
        {
           // If i is more counterclockwise than current q, then
           // update q
           if (orientation(points[p], points[i], points[q]) == 2)
               q = i;
        }
        // Now q is the most counterclockwise with respect to p
        // Set p as q for next iteration, so that q is added to
        // result 'hull'
        p = q;
    } while (p != l);  // While we don't come to first point
    // Print Result
    fout << n << ":";
    cout << n << ":";
    for (int i = hull.size()-1; i>=0;  i--){
        fout << " " <<hull[i].x << " "
              << hull[i].y;
        cout << " " <<hull[i].x << " "
              << hull[i].y;
    }
    fout << "\n";
    cout << "\n";
}
// Driver program to test above functions
int main()
{
    //Read vertex coordinatesfrom Smooth2Polytopes.txt
    //Example vertex coordinates:
    //  [[0,0],[1,0],[0,1]]
    //First gets converted to
    //  @[0,0],[1,0],[0,1]*
    //numbers are always after '[' or after ','
    //end of line is given by two ']'
    ifstream fin("Smooth2Polytopesfixed.txt");
    fout.open("Smooth2Polytopesfixed2.txt");
    char s;
    vector<Point> points= {};
    while(fin >> s){
        int first, second;
        if(s == '['){
            fin >> first >> s >> second;
            Point new_point;
            new_point.x = first;
            new_point.y = second;
            points.push_back(new_point);
        }
        else if(s == '*'){
            int n = points.size();
            convexHull(points, n);
            points = {};
        }
    }
    /*string input_string;
    while (getline(fin, input_string)){
        input_string[0] = '@';
        input_string[input_string.size()-1] = '*';
        fout << input_string << endl;
    }*/
    fout.close();
    return 0;
}
