//====================================================================================
//SOFTWARE:	Network of Eelectrically Conductive Nanocomposites (NECN)
//CODE FILE:	Geometry_2D.cpp
//OBJECTIVE:	The definitions of point, line and shape in 2D
//AUTHOR:		Fei Han
//E-MAIL:			fei.han@kaust.edu.sa
//====================================================================================

#include "Geometry_2D.h"

//---------------------------------------------------------------------------
//The member function for the 2D point class
Point_2D::Point_2D( double px, double py )
{
	x = px;
	y = py;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator+( Point_2D &pt )
{
	Point_2D rp( x + pt.x, y + pt.y);
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator+( const Point_2D &pt )
{
	Point_2D rp( x + pt.x, y + pt.y);
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator+( double d )
{
	Point_2D rp( x + d, y + d);
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator-( Point_2D &pt )
{
	Point_2D rp( x - pt.x, y - pt.y );
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator-( const Point_2D &pt )
{
	Point_2D rp( x - pt.x, y - pt.y );
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator-( double d )
{
	Point_2D rp( x - d, y - d );
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator*( double d )
{
	Point_2D rp( x*d, y*d );
	return rp;
}
//---------------------------------------------------------------------------
Point_2D Point_2D::operator/( double d )
{
	Point_2D rp( x/d, y/d );
	return rp;
}
//---------------------------------------------------------------------------
bool Point_2D::operator==( Point_2D &pt )
{
	return (x==pt.x&&y==pt.y);
}
//---------------------------------------------------------------------------
bool Point_2D::operator!=( Point_2D &pt )
{
	return (x!=pt.x||y!=pt.y);
}
//---------------------------------------------------------------------------
double Point_2D::distance_to(const Point_2D &pt )const
{
	double rv2 = (x-pt.x)*(x-pt.x)+(y-pt.y)*(y-pt.y);
	return sqrt(rv2);
}
//---------------------------------------------------------------------------
double Point_2D::distance_to(const double &px, const double &py)const
{
	double rv2 = (x-px)*(x-px)+(y-py)*(y-py);
	return sqrt(rv2);
}
//===========================================================================

//The member function for the 3D line class
//---------------------------------------------------------------------------
//Constructor
Line_2D::Line_2D(Point_2D p0, Point_2D p1)
{
	point[0] = p0;
	point[1] = p1;
	A = p1.y - p0.y ;
	B = -(p1.x - p0.x);
	D = -(A*p0.x+B*p0.y);
	len = length();
	if(len==0) virtual_line = false;
	else virtual_line = true;
}
//---------------------------------------------------------------------------   
//The length of segment
double Line_2D::length() 
{
	double dx = point[1].x-point[0].x;
	double dy = point[1].y-point[0].y;
	return sqrt(dx*dx+dy*dy);
}
//---------------------------------------------------------------------------
//The distance from a point to a line
double Line_2D::distance_point_to_line(const Point_2D *point_temp)const
{
	if(A==0&&B==0) hout << "Attention: this segment is reduced to a point!" <<endl;
	return fabs(A*point_temp->x+B*point_temp->y+D)/sqrt(A*A+B*B);
}
double Line_2D::distance_point_to_line(const Point_2D &point_temp)const
{
	if(A==0&&B==0) hout << "Attention: this segment is reduced to a point!" <<endl;
	return fabs(A*point_temp.x+B*point_temp.y+D)/sqrt(A*A+B*B);
}
double Line_2D::distance_point_to_line(double dx, double dy)const
{
	if(A==0&&B==0) hout << "Attention: this segment is reduced to a point!" <<endl;
	return fabs(A*dx+B*dy+D)/sqrt(A*A+B*B);
}
//---------------------------------------------------------------------------
//The path from a point to a line, that is, a point is at the side of the line or another side
double Line_2D::path_point_to_line(const Point_2D &point_temp)const
{
	if(A==0&&B==0) hout << "Attention: this segment is reduced to a point!" <<endl;
	return A*point_temp.x+B*point_temp.y+D;
}
//---------------------------------------------------------------------------
//To judge if a point is in a segment
int Line_2D::contain(const Point_2D &point_temp)const
{
	//to judge if the distance from a point to two endpoints is larger than the distance between endpoints
	if( fabs(point_temp.distance_to(point[0])+point_temp.distance_to(point[1])-point[0].distance_to(point[1]))>Zero ) return 0; 
	return 1;
}
int Line_2D::contain(const double &px, const double &py)const
{
	Point_2D point_temp(px, py);
	//to judge if the distance from a point to two endpoints is larger than the distance between endpoints
	if( fabs(point_temp.distance_to(point[0])+point_temp.distance_to(point[1])-point[0].distance_to(point[1]))>Zero ) return 0; 
	return 1;
}
//===========================================================================

//The definition for a rectangle class in 2D
//---------------------------------------------------------------------------
//Constructor
Rectangle::Rectangle(Point_2D p0, Point_2D p1, Point_2D p2, Point_2D p3)
{
	point[0] = p0;
	point[1] = p1;
	point[2] = p2;
	point[3] = p3;
	length = p2.x-p0.x;
	width = p2.y-p0.y;
	area = length*width;
	if(length>0&&width>0) virtual_rect = 1;
	else virtual_rect = 0;
}
Rectangle::Rectangle(Point_2D p0, double len, double wid)
{
	point[0].x = p0.x;		point[0].y = p0.y;
	point[1].x = p0.x;		point[1].y = p0.y+wid;
	point[2].x = p0.x+len;		point[2].y = p0.y+wid;
	point[3].x = p0.x+len;		point[3].y = p0.y;
	length = len;
	width = wid;
	area = length*width;
	if(length>0&&width>0) virtual_rect = 1;
	else virtual_rect = 0;
}
Rectangle::Rectangle(Point_2D p0, Point_2D p2)
{
	point[0] = p0;
	point[1].x = p0.x;
	point[1].y = p2.y;
	point[2] = p2;
	point[3].x = p2.x;
	point[3].y = p0.y;
	length = p2.x-p0.x;
	width = p2.y-p0.y;
	area = length*width;
	if(length>0&&width>0) virtual_rect = 1;
	else virtual_rect = 0;
}
//---------------------------------------------------------------------------
//To calculate the area of the rectangle
double Rectangle::calculate_area()
{
	return length*width;
}
//---------------------------------------------------------------------------
//To judge if a point is in the rectangle (0: no, 1: yes)
int Rectangle::contain(const Point_2D &poi)const
{
	if(  poi.x>=point[0].x && poi.x<=point[2].x &&
		  poi.y>=point[0].y && poi.y<=point[2].y ) return 1;
	else return 0;
}
//---------------------------------------------------------------------------
//To judge if a point is inside the rectanagle but isn't on the boundaries of the rectangle (0: no, 1: yes)
int Rectangle::contain_in(const Point_2D &poi)const
{
	if(  poi.x>point[0].x && poi.x<point[2].x &&
		  poi.y>point[0].y && poi.y<point[2].y ) return 1;
	else return 0;
}
//---------------------------------------------------------------------------
//To judge if a segment is in the rectangle (0: no, 1: yes)
int Rectangle::contain(const Line_2D &line)const
{
	for(int i=0; i<2; i++)
	if( contain(line.point[i])==0 ) return 0;
	return 1;
}
//---------------------------------------------------------------------------
//To judge if a rectangle is in this rectangle (0: no, 1: yes)
int Rectangle::contain(const Rectangle &rect)const
{
	for(int i=0; i<4; i++)
	if( contain(rect.point[i])==0 ) return 0;
	return 1;
}
//---------------------------------------------------------------------------
//To judge if a segment is inside the rectangle, but doesn't intersect with the boundaries of the rectangle (0: no, 1: yes)
int Rectangle::contain_in(const Line_2D &line)const
{
	if( contain(line.point[0])==1&&contain(line.point[1])==1 ) 
	{
		if( line.point[0].x == line.point[1].x&&(line.point[0].x==point[0].x||line.point[0].x==point[2].x) ) return 0;
		if( line.point[0].y == line.point[1].y&&(line.point[0].y==point[0].y||line.point[0].y==point[2].y) ) return 0;
		return 1;
	}
	return 0;
}
//---------------------------------------------------------------------------
//To judge if a segment is intersect with the rectangle, including both endpoints of the segments on the boundaries of the rectangle or at least one endpoint inside the rectangle.
//If only one endpoint at the boundaries of the rectangle but another endpoint is outside of the rectangle, it isn't considered as intersection) (0: no, 1: yes)
int Rectangle::overlap(const Line_2D &line)const
{
	if( contain(line.point[0])==1||contain(line.point[1])==1 )
	{
		if( line.point[0].x == line.point[1].x) 
		{
			if( line.point[0].y > line.point[1].y&&(line.point[0].y <= point[0].y||line.point[1].y>=point[2].y) ) return 0;
			if( line.point[0].y < line.point[1].y&&(line.point[1].y <= point[0].y||line.point[0].y>=point[2].y) ) return 0;
		}
		if( line.point[0].y == line.point[1].y) 
		{
			if( line.point[0].x > line.point[1].x&&(line.point[0].x <= point[0].x||line.point[1].x>=point[2].x) ) return 0;
			if( line.point[0].x < line.point[1].x&&(line.point[1].x <= point[0].x||line.point[0].x>=point[2].x) ) return 0;
		}
		return 1;
	}
	return 0;
}
//---------------------------------------------------------------------------
//Calculate the overlapping rectangle between Rect1 and Rect2
int Rectangle::rectangles_overlapping(const Rectangle &Rect1, const Rectangle &Rect2)
{
	//the common part of both rectangles, that is, to sort the points positions from lower to upper, and then to choose the middle two points
	//to choose the point with the bigger value in both lower points of two rectangles
	if(Rect1.point[0].x<=Rect2.point[0].x)	point[0].x = Rect2.point[0].x;
	else point[0].x = Rect1.point[0].x;
	if(Rect1.point[0].y<=Rect2.point[0].y)	point[0].y = Rect2.point[0].y;
	else point[0].y = Rect1.point[0].y;
	//to choose the point with the smaller value in both upper points of two rectangles
	if(Rect1.point[2].x<=Rect2.point[2].x)	point[2].x = Rect1.point[2].x;
	else point[2].x = Rect2.point[2].x;
	if(Rect1.point[2].y<=Rect2.point[2].y)	point[2].y = Rect1.point[2].y;
	else point[2].y = Rect2.point[2].y;
	
	point[1].x = point[0].x;
	point[1].y = point[2].y;
	point[3].x = point[2].x;
	point[3].y = point[0].y;
	length = point[2].x-point[0].x;
	width = point[2].y-point[0].y;
	area = length*width;
	if(length>0&&width>0) virtual_rect = 1;
	else virtual_rect = 0;  //Two rectangles have no overlapping, or only overlap with a vertex or a boundary

	return virtual_rect;
}
//---------------------------------------------------------------------------
//Divide the rectangle into a Sudoku where the inner rectangle becomes the center element of the Sudoku (将矩形以内部小矩形为中心做九宫格分割)
int Rectangle::make_nine_grids_by(const Rectangle rect, vector<Rectangle> &grids)
{
	if(contain(rect)==0) 
	{
		hout << "This rectangle is not inside the rectangle which will be divided, please check it!" << endl;
		return 0;
	}	
	//to record the basic points for division
	Point_2D point_base[4];
	point_base[0] = point[0];
	point_base[1] = rect.point[0];
	point_base[2] = rect.point[2];
	point_base[3] = point[2];
	//to calculate all points on the Sudoku (from left to right, from lower to upper)
	Point_2D point_ex[16];
	for(int i=0; i<4; i++)
		for(int j=0; j<4; j++)
		{
			point_ex[i*4+j].x = point_base[j].x;
			point_ex[i*4+j].y = point_base[i].y;
		}
	//to record every sub-rectangle in the Sudoku
	grids.clear();
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
		{
			Rectangle rect_temp( point_ex[i*4+j], point_ex[(i+1)*4+j], point_ex[(i+1)*4+j+1], point_ex[i*4+j+1] );
			grids.push_back(rect_temp);
		}

	return 1;
}

//---------------------------------------------------------------------------
//Output all information of a rectangle, including four vertex, length, width, area and virtual_rect
void Rectangle::output_parameters()
{
	hout << "The coordinates of four corner points (starting from lower-left corner point in clockwise direction):" << endl;
	for(int i=0; i<4; i++) hout << "(" << point[i].x << ", " << point[i].y << ")" << "   ";
	hout << endl;
	hout << "The length, width, area and virtual_rect tag:" << endl;
	hout << "Length " << length << ", Width " << width << ", area " << area << ", virtual_rect " << virtual_rect << endl;
}
//===========================================================================
