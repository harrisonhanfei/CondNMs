//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Geometry_2D.h
//OBJECTIVE:	The definitions of point, line and shape in 2D
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#ifndef GEOMETRY_2D_H
#define GEOMETRY_2D_H

#include<cmath>
#include<stdlib.h>
#include<vector>
#include "Hns.h"
using namespace hns;

#include "MathMatrix.h"

//---------------------------------------------------------------------------
//The definition for points in 2D
class Point_2D
{
	public:

		//Data Member
		double x, y;
		int flag;

		//Constructor
		Point_2D(){};
		Point_2D( double px, double py );

		//Member Functions
        Point_2D operator+( Point_2D &pt );
        Point_2D operator+( const Point_2D &pt );
        Point_2D operator+( double d );
        Point_2D operator-( double d );
        Point_2D operator-( Point_2D &pt );
        Point_2D operator-( const Point_2D &pt );
        Point_2D operator*( double d );
        Point_2D operator/( double d );
		bool operator==( Point_2D &pt );
		bool operator!=( Point_2D &pt );
        double distance_to(const Point_2D &pt)const;
		double distance_to(const double &px, const double &py)const;
};

//---------------------------------------------------------------------------
//The definition for lines (segements) in 2D
class Line_2D
{
	public:

		//Data Member
		 Point_2D point[2];		//the coordinates of two endpoints of a segment
		 double len;					//the length of a segment
		 bool virtual_line;			//to mark if it is a virtual(false) segment (false: reduced to a point; true: a real segment)


		//Constructor
		Line_2D(){};
		Line_2D(Point_2D p0, Point_2D p1);
		
		//Member Functions
		double length();																						//the length of a segment
        double distance_point_to_line(const Point_2D *point_temp)const;			//the distance from a point to a line
        double distance_point_to_line(const Point_2D &point_temp)const;		//the distance from a point to a line
        double distance_point_to_line(double dx, double dy)const;					//the distance from a point to a line
        double path_point_to_line(const Point_2D &point_temp)const;				//the path from a point to a line, that is, a point is at the side of the line or another side
		int contain(const Point_2D &point_temp)const;										//to judge if a point is in a segment
		int contain(const double &px, const double &py)const;							//to judge if a point is in a segment

	private:
				
		//Data Member
        double A, B, D ;       //the coeffecients of a line equation: Ax+By+D=0
};

//---------------------------------------------------------------------------
//The definition for a rectangle class in 2D
class Rectangle
{
	public:

		//Data Member
		Point_2D point[4];		//the four vertex of the rectangle, account from left bottom in a clockwise direction
		double length;				//the length of the rectangle
		double width;				//the width of the rectangle
		double area;					//the area of the rectangle
		int virtual_rect;				//to mark if it is a virtual(false) rectangle (false: not a rectangle; true: a real rectangle)

		//Constructor
		Rectangle(){};
		Rectangle (Point_2D p0, Point_2D p1, Point_2D p2, Point_2D p3);
		Rectangle (Point_2D p0, double len, double wid);
		Rectangle (Point_2D p0, Point_2D p2);

		//Member Functions
		double calculate_area();									//to calculate the area of the rectangle
		int contain(const Point_2D &poi)const;			//to judge if a point is in the rectangle (0: no, 1: yes)
		int contain_in(const Point_2D &poi)const;		//to judge if a point is inside the rectanagle but isn't on the boundaries of the rectangle (0: no, 1: yes)
		int contain(const Line_2D &line)const;			//to judge if a segment is in the rectangle (0: no, 1: yes)
		int contain(const Rectangle &rect)const;		//to judge if a rectangle is in this rectangle (0: no, 1: yes)
		int contain_in(const Line_2D &line)const;		//to judge if a segment is inside the rectangle, but doesn't intersect with the boundaries of the rectangle (0: no, 1: yes)
		int overlap(const Line_2D &line)const;			//to judge if a segment is intersect with the rectangle, including both endpoints of the segments on the boundaries of the rectangle 
																			//or at least one endpoint inside the rectangle.
																			//If only one endpoint at the boundaries of the rectangle but another endpoint is outside of the rectangle, it isn't considered as intersection) (0: no, 1: yes)
		void output_parameters();								//Output all information of a rectangle, including four vertex, length, width, area and virtual_rect
		int rectangles_overlapping(const Rectangle &Rect1, const Rectangle &Rect2);			//Calculate the overlapping rectangle between Rect1 and Rect2
		int make_nine_grids_by(const Rectangle rect, vector<Rectangle> &grids);					//Divide the rectangle into a Sudoku where the inner rectangle becomes the center element of the Sudoku (将矩形以内部小矩形为中心做九宫格分割)
};
//---------------------------------------------------------------------------
#endif
//===========================================================================
