//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Contact_grid.cpp
//OBJECTIVE:	Create a background grid to find the points in contact faster
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#include "Contact_grid.h"

//	This function groups the points inside the observation window into sub-regions.
//	The task of finding contacting points can be computationally expensive. To reduce computational cost, the sample was divided into smaller cubic sub-regions. 
//	The cubic sample of size $[a \times a \times a]$ was divided into $m$ segments in each direction resulting in $m^3$ cubic sub-regions. 
//	Then, the contacting points are searched for only on each smaller sub-region rather than in the whole sample. It may happen that two contacting points belong to different sub-regions. 
//	In order to take into account these contacting points, each sub-region is ''extended". If a point lies exactly at a boundary, then a point can be at a maximum distance equal to
//	$r_{max}+r_{max}+d_W = 2r_{max}+d_W$ from the boundary. 
//	Here, $r_{max}$ is the maximum radii of the CNTs. Then, each boundary plane of each cubic sub-region is translated externally to a parallel plane at a distance $2r_{max}+d_W$.
 
// Input:
//    struct Geom_RVE sample
//        Goemetry of the simulated sample
//    struct Cutoff_dist cutoffs
//        Structure that contains the cutoff for tunneling and overlapping
//    struct Nanotube_Geo cnts
//        Structure that contains the geometry of the CNTs
//    vector<int> cnts_inside
//        List of CNTs that are inside the observation window
//    vector<vector<long int> > structure
//        Vector with the structure
//   vector<Point_3D> points_in
//        List of points
//    int window
//        Current observation window. window=0 means the largest observation window
 
// Output:
//    vector<vector< long int> > sectioned_domain
//       Vector with overlaping sub-regions. This vector is used to look for contact points faster
 
int Contact_grid::Generate_contact_grid(const struct Geom_RVE &sample, const struct Cutoff_dist &cutoffs, const struct Nanowire_Geo &nanowires, const vector<int> &cnts_inside, 
															 const vector<vector<long int> > &structure, vector<Point_3D> &points_in, int window)
{
    //Dimensions of the current observation window
    double w_x = sample.win_max_x - window*sample.win_delt_x;
    double w_y = sample.win_max_y - window*sample.win_delt_y;
    double w_z = sample.hei_z;

    //These variables are the coordinates of the lower corner of the observation window
    double xmin = sample.origin.x + (sample.len_x - w_x)/2;
    double ymin = sample.origin.y + (sample.wid_y - w_y)/2;
    double 	zmin = sample.origin.z;
    
    //Sizes of each region
    double dx = sample.gs_minx;
    double dy = sample.gs_miny;
    double dz = sample.hei_z;
    
    //Number of regions on each direction
    int sx = (int)(w_x/dx);
    int sy = (int)(w_y/dy);

    //Maximum distance between two points in contact inside the sample
    double cutoff = cutoffs.tunneling_dist + 2*nanowires.rad_max;

    //Check that the regions are not too small for the maximum cutoff distance 2r_max+tunnel
    //If they are, then change the number of sections to the maximum possible
    if (dx < 2*cutoff) 
	{
        sx = (int)(w_x/(2*(cutoff+Zero)));
        dx = w_x/(double)sx;
        hout << "Modified the number of sections along x. " << "sx = " << sx << '\t' << "dx = " << dx << endl;
    }
    if (dy < 2*cutoff)
	{
        sy = (int)(w_y/(2*(cutoff+Zero)));
        dy = w_y/(double)sy;
        hout << "Modified the number of sections along y. " << "sy = " << sy << '\t' << "dy = " << dy << endl;
    }

	//These variables will give me the region cordinates of the region that a point belongs to
    int a, b;
    int t;
    
    //These variables are to reduce operations when accessing the coordinates of a point and it's CNT number
    double x, y;
    int CNT;
    long int P;
    
    //There will be sx*sy different regions
    sectioned_domain.clear();
    vector<long int> empty_long;
    sectioned_domain.assign(sx*sy, empty_long);
    hout << "There are " << sx*sy << " overlapping sub-regions." << endl;
    
    //First loop over the CNTs inside the box, then loop over the points inside each CNT
    for (int i = 0; i < (int)cnts_inside.size(); i++) 
	{
        CNT = cnts_inside[i];
        for (int j = 0; j < (int)structure[CNT].size(); j++) 
		{
            P = structure[CNT][j];
            //Save coordinates of the point
            x = points_in[P].x;
            y = points_in[P].y;
            
            //Calculate the region-coordinates
            a = (int)((x-xmin)/dx);
            //Limit the value of a as it has to go from 0 to sx-1
            if (a == sx) a--;

            b = (int)((y-ymin)/dy);
            //Limit the value of b as it has to go from 0 to sy-1
            if (b == sy) b--;

            //Coordinates of non-overlaping region the point belongs to
            double x1 = a*dx +  xmin;
            double x2 = x1 + dx;
            double y1 = b*dy +  ymin;
            double y2 = y1 + dy;
            
            //Initialize flags for overlaping regions
            int fx = 0;
            int fy = 0;
            
            //Assign value of flag according to position of point
            //The first operand eliminates the periodicity on the boundary
            if ((x > cutoff + xmin) && (x >= x1) && (x <= x1+cutoff))
                fx = -1;
            else if ((x < w_x+xmin-cutoff) && (x >= x2-cutoff) && (x <= x2 ))
                fx = 1;
            if ((y > cutoff + ymin) && (y >= y1) && (y <= y1+cutoff))
                fy = -1;
            else if ((y < w_y+ymin-cutoff) && (y >= y2-cutoff) && (y <= y2 ))
                fy = 1;
            
            //Create array for loop over overlaping regions
            int temp[2][2] = {{a+fx, b+fy}, {a, b}};
            
            //In this loop I check all regions a point can belong to when it is in an overlaping zone
            for (int ii = 0; ii < 2; ii++) 
			{
                if (!fx) ii++;							//if flag is zero, do this loop only once
                for (int jj = 0; jj < 2; jj++) 
				{
                    if (!fy) jj++;						//if flag is zero, do this loop only once
					t = calculate_t_2D(temp[ii][0], temp[jj][1], sx);
					sectioned_domain[t].push_back(P);
                }
            }
        }
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//Calculates the region to which a point corresponds
int Contact_grid::calculate_t(int a, int b, int c, int sx, int sy)
{
    return a + b*sx + c*sx*sy;
}
//---------------------------------------------------------------------------
//Calculates the 2D region to which a point corresponds
int Contact_grid::calculate_t_2D(int a, int b, int sx)
{
    return a + b*sx;
}
//===========================================================================
