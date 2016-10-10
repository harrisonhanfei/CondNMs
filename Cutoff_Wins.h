//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Cutoff_Wins.h
//OBJECTIVE:	To cutoff the windows
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef CUTOFFWINS_H
#define CUTOFFWINS_H

#include "Input_Reader.h"
#include "Background_vectors.h"
#include "Tecplot_Export.h"

//-------------------------------------------------------
class Cutoff_Wins
{
public:
    //Data Member
    vector<int> cnts_inside;								//List of CNTs inside the observation window
    vector<vector<int> > boundary_cnt;				//Boundary vector. It is used to determine percolation
	vector<vector<short int> > boundary_flags;	//This vector will help find points on the boundary. This is used in the direct electrifying algorithm

    //Constructor
    Cutoff_Wins(){};

    //Member Functions
    int Extract_observation_window(const struct Geom_RVE sample, const struct Nanowire_Geo cnts, vector<vector<long int> > &structure, vector<double> &radii, vector<Point_3D> &points_in, vector<vector<int> > &shells_cnt, const int &window);
    
private:
    vector<vector<int> > sectioned_domain;
    double xmin, ymin, zmin;
    double w_x, w_y, w_z;    
	
	//Constructor

	//Member Functions
	//To scan every nanotube that is in the boundary region and to delete and to trim CNTs outside.
    int Trim_boundary_cnts(vector<vector<int> > &shells_cnt, const int &window, const struct Geom_RVE &sample, vector<Point_3D> &points_in, vector<vector<long int> > &structure, vector<double> &radii);
    int First_index(vector<Point_3D> &points_in, const vector<long int> &structure_CNT, const int &new_CNT, int &index1);
    int Second_index(vector<Point_3D> &points_in, vector<long int> &structure_CNT, int &new_CNT, int &index2);
	//(Quasi2D Ag nanowire) this function checks the three location in which a point is placed
    string Where_is(const Point_3D &point)const;
    int New_boundary_point(struct Geom_RVE sample, vector<Point_3D> &points_in, vector<vector<long int> > &structure, long int insidePoint, long int outsidePoint, int CNT, string currentLocation);
    int Substitute_boundary_point(vector<Point_3D> &points_in, long int global_i, long int global_o);
    int Get_intersecting_point_RVE_surface_quasi2D(const Point_3D &point0, const Point_3D &point1, vector<Point_3D> &ipoi_vec);
	//(Quasi2D)Add the corrent point to the corrsponding boundary vector
    void Add_to_boundary_vectors(const Point_3D &point3d, const long int &point, const int &new_CNT);
	//This function adds a CNT to the corresponding boundary vector
    void Add_CNT_to_boundary(vector<int> &boundary, const int &CNT, const long int &point, const short int &flag1, const short int &flag2);    

};
//-------------------------------------------------------
#endif
//===========================================================================