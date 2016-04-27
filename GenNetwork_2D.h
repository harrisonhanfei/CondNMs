//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	GenNetwork.h
//OBJECTIVE:	To generate networks with overlapping
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef GENNETWORK2D_H
#define GENNETWORK2D_H

#include<iomanip>
#include<iostream>
#include<sstream>
#include<fstream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<math.h>
#include <random> 
#include "Hns.h"
#include "Input_Reader.h"
#include "Geometry_3D.h"
#include "MathMatrix.h"
#include "Fem_3D.h"
#include "Gauss.h"
#include "Tecplot_Export.h"

using namespace hns;

//-------------------------------------------------------
class GenNetwork
{
	public:
		//Data Member

		//Constructor
		GenNetwork(){};

		//Member Functions
		int Generate_nanowire_networks(const struct Geom_RVE &geom_rve, const struct Nanowire_Geo &nanowire_geo, vector<Point_3D> &cpoints, vector<double> &cnts_radius, vector<vector<long int> > &cstructures)const;
		//Generate the nodes and tetrahedron elements of nanowires (No const following this function because a sum operation on two Point_2D points inside)
		int Generate_cnts_nodes_elements(vector<vector<Node> > &nodes, vector<vector<Element> > &eles, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radius);
        //Generate the nodes and tetrahedron elements of nanotubes (No const following this function because a sum operation on two Point_3D points inside). This function uses a 1D point vector and a 2D structure vector that references the point vector
        int Generate_cnts_nodes_elements(vector<vector<Node> > &nodes, vector<vector<Element> > &eles, const vector<Point_3D> &cnts_points, const vector<double> &cnts_radius, const vector<vector<long int> > &structure);
	private:
		//Data Member

		//Member Functions

		//Generate a network defined by points and connections 
		int Generate_network_threads_mt(const struct Geom_RVE &geom_rve, const struct Nanowire_Geo &nanowire_geo, vector<vector<Point_3D> > &Cnts_points, vector<double> &Cnts_radius)const;
		//---------------------------------------------------------------------------
		//---------------------------------------------------------------------------
		//---------------------------------------------------------------------------
		int Get_random_value_mt(const string &dist_type, mt19937 &engine, uniform_real_distribution<double> &dist, const double &min, const double &max, double &value)const;
		int Get_seed_point_mt(const struct cuboid &cub, Point_3D &point, mt19937 &engine_x, mt19937 &engine_y, uniform_real_distribution<double> &dist)const;
		int Get_uniform_direction_mt(const struct Nanowire_Geo &nanowire_geo, double &cnt_pha, mt19937 &engine_pha, uniform_real_distribution<double> &dist)const;
		//---------------------------------------------------------------------------
		//---------------------------------------------------------------------------
		//---------------------------------------------------------------------------
		//To judge if a point is included in a RVE
		int Judge_RVE_including_point(const struct cuboid &cub, const Point_3D &point)const;
        //To judge if a point is included in a RVE
        int Judge_RVE_including_point(const struct Geom_RVE &geom_rve, const Point_3D &point)const;
		//Calculate all intersection points between the new segment and surfaces of RVE
		//(using a parametric equatio:  the parameter 0<t<1, and sort all intersection points from the smaller t to the greater t)  
		int Get_intersecting_point_RVE_surface(const struct cuboid &cub, const Point_3D &point0, const Point_3D &point1, vector<Point_3D> &ipoi_vec)const;
		//This functions initializes vector sub-regions
		void Initialize_subregions(const struct Geom_RVE &geom_rve, vector<int> &nsubregions, vector<vector<long int> > &sectioned_domain)const;
		//---------------------------------------------------------------------------
		//---------------------------------------------------------------------------
		//---------------------------------------------------------------------------
		//To calculate the effective portion (length) which falls into the given region (RVE)
		int Effective_length_given_region(const struct cuboid &cub, const Point_3D &last_point, const Point_3D &new_point, double &real_length)const;
		//Transform the 2D cnts_points into 1D cpoints and 2D cstructuers
		int Transform_cnts_points(const vector<vector<Point_3D> > &cnts_points, vector<Point_3D> &cpoints, vector<vector<long int> > &cstructures)const;
		//---------------------------------------------------------------------------
		//---------------------------------------------------------------------------
		//---------------------------------------------------------------------------
		//Transform angles into matrix
		MathMatrix Get_transformation_matrix(const double &sita, const double &pha)const;
		int Get_angles_vector_in_spherial_coordinates(const Point_3D &normal, double &sita, double &pha)const;
		//Calculate a group of equidistant points along the circumference which is on the plane defined by the center point of the circle and the normal vector
		int Get_points_circle_in_plane(const Point_3D &center, const double &sita, const double &pha, const double &radius, const int &num_sec, vector<Node> &nod_temp)const;
		//Calculate a group of projected points (which are on the plane with the center point of the circle and the normal vector) 
		//which are projected from a group of points on the previous circumference and projected along the direction of line_vec
		int Get_projected_points_in_plane(const Point_3D &center, const Point_3D &normal, const Point_3D &line, const int &num_sec, vector<Node> &nod_temp)const;
};
//-------------------------------------------------------
#endif
//===========================================================================