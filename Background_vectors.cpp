//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Background_vectors.cpp
//OBJECTIVE:	Using nested shells on the background to mark the CNTs for trimming faster in each observation window
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Background_vectors.h"


 // This class generates the vector:
 //	shells_cnt: is used to trim the CNTs. The CNTs are grouped according to the nested regions, then need to be trimmed
 //	The shells will be created as follows: the smallest observation window (cuboid) will be a "core" region (smallest shell) and will be used for the last element in the vector.
 //	Then, the next shell-subregion will be the volume of the next observation window minus the volume of the first one (the core region mentioned above).
 //	The following shells will have the same form: the volume of the observation window minus the volume of the previous one.
 //	The last shell region will be the boundary layer. This will be used for the first element of the vector
 //
 //	Input:
 //   struct Geom_RVE sample
 //       Geometry of the generated sample
 //   struct Nanotube_Geo cnts
 //       Geometry of the CNTs
 //   vector<Point_3D> points_in
 //       List of points
 //
 //	Output:
 //   vector<vector<int> > shells_cnt
 //
 //	Modified inputs:

//---------------------------------------------------------------------------
//Generate the shell vector
int Background_vectors::Generate_shells_and_structure(const struct Geom_RVE &sample, const struct Nanowire_Geo &cnts, const vector<Point_3D> &points_in, vector<vector<int> > &shells_cnt)const
{
    //Initialize the shells_cnt vector
    //Empty vector to initialize the shells_cnt vector
    vector<int> empty;
    //sample.cut_num is the number of increments from the smallest to the largest observation widows
    //Then, there are sample.cut_num+1 observation windows
    //Then, there are sample.cut_num+2 shell sub-regions if we include the boundary layer
    //Hence the shells_cnt vector will have sample.cut_num+2 elements
    shells_cnt.assign(sample.cut_num+2, empty);
    //hout << "sample.cut_num+2="<<sample.cut_num+2<<endl;
    
    //Scan all points to determine in which sub-regions the CNTs are located and construct the structure vector
    for(long int i = 0; i<(long int)points_in.size(); i++)
	{
        //Add to the corresponding shell
        if(!Add_to_shell(sample, points_in[i], shells_cnt)) 
		{
            hout << "Error happens at Add_to_shell() in the Generate_shells_and_structure() function."<< endl;
            return 0;
        }
    }

    return 1;
}

//---------------------------------------------------------------------------
//This function finds the corresponding shell sub-region where the point is located. Then it adds the CNT number to the
//corresponding vector in the 2D vector shells_cnt (for the quasi-2D network of Ag nanowires )
int Background_vectors::Add_to_shell(const struct Geom_RVE &sample, const Point_3D &point, vector<vector<int> > &shells_cnt)const
{
    //Find the shell based on the x coordinate
    int shell_x = Find_shell(point.x, sample.origin.x, sample.len_x, sample.win_delt_x, sample.win_min_x, sample.win_max_x, shells_cnt);
    //Find the shell based on the y coordinate
    int shell_y = Find_shell(point.y, sample.origin.y, sample.wid_y, sample.win_delt_y, sample.win_min_y, sample.win_max_y, shells_cnt);
	if(shell_x==-1||shell_y==-1) 
	{
		hout << "Error: shell_x=" << shell_x << " shell_y=" << shell_y << " calculated by Find_shell() in the Add_to_shell() function."<< endl;
		return 0;
	}

    //The shell to which the CNT of point is the outer-most shell from the three coordinates, that is, the minimum shell number overall
    int shell = shell_x;
    //With this if-statement I have shell = min(shell_x, shell_y)
    if(shell_y<shell) shell = shell_y;
            
    //Finally add the CNT on the corresponding sectioned domain
    //Add the CNT number to the shell sub-region, only if it is empty or if the last CNT is not the current CNT
    if((!shells_cnt[shell].size()) || (shells_cnt[shell].back() != point.flag)) shells_cnt[shell].push_back(point.flag);
    
    return 1;
}

//---------------------------------------------------------------------------
//This function finds the shell to which one coordinate belongs to
int Background_vectors::Find_shell(const double &x_in, const double &x_0, const double &len_x, const double &dx, const double &win_min_x, const double &win_max_x, vector<vector<int> > &shells_cnt)const
{
    //Calculate the middle point of the sample
    double x_m = x_0 + len_x/2;
    //Calculate the minimum coordinate of the largest observation window
    //below this corrdinate, x is outside the largest observation window,
    //or in the boundary layer if the largest observation window equals the sample size
    double x_layer = x_0 + (len_x - win_max_x)/2;
    //The minimum coordinate of the smallesr observation window
    //any point above this coordinate is inside the core shell
    double x_core = x_0 + (len_x - win_min_x)/2;
    //An observation window grows a total of dx. However, as it is centered, it gros dx/2 on each direction
    double hfdx = dx/2;
    
    //Effective coordinate
    double x;
    //If the point is greater than the middle point, map it to the mirrored range
    if (x_in > x_m) 
	{
        double m = (x_m - x_0)/(x_m - (x_0 + len_x)); //x_0 + len_x = maximum coordinate of the sample
        x = m*(x_in - x_m) + x_m;
    } 
	else x = x_in;

    //Check if x is in the outer shell (it is in the outer shell when the point is below the x_min)
    if (x < x_layer)  return 0;
	else if( x > x_core )  return ((int)shells_cnt.size()-1);	//Check if x is in the inner shell
	else return ((int)((x - x_layer + Zero)/hfdx)+1);			//The Zero is for foating point errors
    
    return -1;	//Error happens when reaches here
}

//===========================================================================
