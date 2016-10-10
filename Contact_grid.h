//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Contact_grid.h
//OBJECTIVE:	Create a background grid for faster finding the points in contact
//AUTHOR:		Fei Han; Angel Mora
//E-MAIL:			fei.han@kaust.edu.sa	;	angel.mora@kaust.edu.sa
//====================================================================================

#ifndef CONTACT_GRID_H
#define CONTACT_GRID_H

#include <iostream>
#include "Input_Reader.h"

//-------------------------------------------------------
class Contact_grid
{
public:
    //Data Member
    vector<vector< long int> > sectioned_domain;
    
    //Constructor
    Contact_grid(){};
    
    //Member Functions
    int Generate_contact_grid(const struct Geom_RVE &sample, const struct Cutoff_dist &cutoffs, const struct Nanowire_Geo &nanowires,
											const vector<int> &cnts_inside, const vector<vector<long int> > &structure, vector<Point_3D> &points_in, int window);
    int calculate_t(int a, int b, int c, int sx, int sy);
    int calculate_t_2D(int a, int b, int sx);

private:

};
//-------------------------------------------------------
#endif
//===========================================================================
