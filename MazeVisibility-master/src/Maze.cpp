/************************************************************************
	 File:        Maze.cpp

	 Author:
				  Stephen Chenney, schenney@cs.wisc.edu
	 Modifier
				  Yu-Chi Lai, yu-chi@cs.wisc.edu

	 Comment:
						(c) 2001-2002 Stephen Chenney, University of Wisconsin at Madison

						Class header file for Maze class. Manages the maze.


	 Platform:    Visio Studio.Net 2003 (converted to 2005)

*************************************************************************/

#include "Maze.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <FL/Fl.h>
#include <FL/fl_draw.h>
#include <GL/GL.h>
#include <GL/glu.h>
#include<windows.h>
const char Maze::X = 0;
const char Maze::Y = 1;
const char Maze::Z = 2;

const float Maze::BUFFER = 0.1f;
const bool Display = false;

void SetColor(int color = 7)
{
	HANDLE hConsole;
	hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(hConsole, color);
}

//**********************************************************************
//
// * Constructor for the maze exception
//======================================================================
MazeException::
MazeException(const char* m)
//======================================================================
{
	message = new char[strlen(m) + 4];
	strcpy(message, m);
}


//**********************************************************************
//
// * Constructor to create the default maze
//======================================================================
Maze::
Maze(const int nx, const int ny, const float sx, const float sy)
//======================================================================
{
	// Build the connectivity structure.
	Build_Connectivity(nx, ny, sx, sy);

	// Make edges transparent to create a maze.
	Build_Maze();

	// Set the extents of the maze
	Set_Extents();

	// Default values for the viewer.
	viewer_posn[X] = viewer_posn[Y] = viewer_posn[Z] = 0.0;
	viewer_dir = 0.0;
	viewer_fov = 45.0;

	// Always start on the 0th frame.
	frame_num = 0;
}


//**********************************************************************
//
// * Construtor to read in precreated maze
//======================================================================
Maze::
Maze(const char* filename)
//======================================================================
{
	char    err_string[128];
	FILE* f;
	int	    i;

	// Open the file
	if (!(f = fopen(filename, "r")))
		throw new MazeException("Maze: Couldn't open file");

	// Get the total number of vertices
	if (fscanf(f, "%d", &num_vertices) != 1)
		throw new MazeException("Maze: Couldn't read number of vertices");

	// Read in each vertices
	vertices = new Vertex * [num_vertices];
	for (i = 0; i < num_vertices; i++) {
		float x, y;
		if (fscanf(f, "%g %g", &x, &y) != 2) {
			sprintf(err_string, "Maze: Couldn't read vertex number %d", i);
			throw new MazeException(err_string);
		}
		vertices[i] = new Vertex(i, x, y);
	}

	// Get the number of edges
	if (fscanf(f, "%d", &num_edges) != 1)
		throw new MazeException("Maze: Couldn't read number of edges");

	// read in all edges
	edges = new Edge * [num_edges];
	int max_horizon = 0;
	int max_vertical = 0;
	for (i = 0; i < num_edges; i++) {
		int     vs, ve, cl, cr, o;
		float	r, g, b;
		if (fscanf(f, "%d %d %d %d %d %g %g %g",
			&vs, &ve, &cl, &cr, &o, &r, &g, &b) != 8) {
			sprintf(err_string, "Maze: Couldn't read edge number %d", i);
			throw new MazeException(err_string);
		}
		edges[i] = new Edge(i, vertices[vs], vertices[ve], r, g, b);
		edges[i]->Add_Cell((Cell*)cl, Edge::LEFT);
		edges[i]->Add_Cell((Cell*)cr, Edge::RIGHT);
		edges[i]->opaque = o ? true : false;

		if (max_horizon < this->edges[i]->endpoints[Edge::START]->posn[Vertex::X])
		{
			max_horizon = this->edges[i]->endpoints[Edge::START]->posn[Vertex::X];
		}
		if (max_horizon < this->edges[i]->endpoints[Edge::END]->posn[Vertex::X])
		{
			max_horizon = this->edges[i]->endpoints[Edge::END]->posn[Vertex::X];
		}
		if (max_vertical < this->edges[i]->endpoints[Edge::START]->posn[Vertex::Y])
		{
			max_vertical = this->edges[i]->endpoints[Edge::START]->posn[Vertex::Y];
		}
		if (max_vertical < this->edges[i]->endpoints[Edge::END]->posn[Vertex::Y])
		{
			max_vertical = this->edges[i]->endpoints[Edge::END]->posn[Vertex::Y];
		}
	}
	//max_horizon -= 2; max_vertical -= 2;
	max_radius = ceil(sqrt(pow(max_horizon, 2) + pow(max_vertical, 2)));
	//printf("%d\n", max_radius);
	my_edge = new double* [num_edges];
	isUsed = new int[num_edges];
	isTransperant = new int[num_edges];
	for (int i = 0; i < (int)this->num_edges; i++)
	{
		my_edge[i] = new double[7];
		my_edge[i][0] = (double)this->edges[i]->endpoints[Edge::START]->posn[Vertex::X];
		my_edge[i][1] = (double)this->edges[i]->endpoints[Edge::START]->posn[Vertex::Y];
		my_edge[i][2] = (double)this->edges[i]->endpoints[Edge::END]->posn[Vertex::X];
		my_edge[i][3] = (double)this->edges[i]->endpoints[Edge::END]->posn[Vertex::Y];
		my_edge[i][4] = (double)this->edges[i]->color[0];
		my_edge[i][5] = (double)this->edges[i]->color[1];
		my_edge[i][6] = (double)this->edges[i]->color[2];
		if (Display)
		{
			printf("%d: ( %f, %f)->( %f, %f),(%f,%f,%f), ", i, my_edge[i][0], my_edge[i][1], my_edge[i][2], my_edge[i][3], 255 * my_edge[i][4], 255 * my_edge[i][5], 255 * my_edge[i][6]);
			if (edges[i]->opaque)
			{
				printf("opaque\n");
			}
			else
			{
				printf("transparent\n");
			}
		}
	}
	// Read in the number of cells
	if (fscanf(f, "%d", &num_cells) != 1)
		throw new MazeException("Maze: Couldn't read number of cells");


	// Read in all cells
	cells = new Cell * [num_cells];
	for (i = 0; i < num_cells; i++) {
		int epx, epy, emx, emy;
		if (fscanf(f, "%d %d %d %d", &epx, &epy, &emx, &emy) != 4) {
			sprintf(err_string, "Maze: Couldn't read cell number %d", i);
			throw new MazeException(err_string);
		}
		cells[i] = new Cell(i, epx >= 0 ? edges[epx] : NULL,
			epy >= 0 ? edges[epy] : NULL,
			emx >= 0 ? edges[emx] : NULL,
			emy >= 0 ? edges[emy] : NULL);
		if (cells[i]->edges[0]) {
			if (cells[i]->edges[0]->neighbors[0] == (Cell*)i)
				cells[i]->edges[0]->neighbors[0] = cells[i];
			else if (cells[i]->edges[0]->neighbors[1] == (Cell*)i)
				cells[i]->edges[0]->neighbors[1] = cells[i];
			else {
				sprintf(err_string,
					"Maze: Cell %d not one of edge %d's neighbors",
					i, cells[i]->edges[0]->index);
				throw new MazeException(err_string);
			}
		}

		if (cells[i]->edges[1]) {
			if (cells[i]->edges[1]->neighbors[0] == (Cell*)i)
				cells[i]->edges[1]->neighbors[0] = cells[i];
			else if (cells[i]->edges[1]->neighbors[1] == (Cell*)i)
				cells[i]->edges[1]->neighbors[1] = cells[i];
			else {
				sprintf(err_string,
					"Maze: Cell %d not one of edge %d's neighbors",
					i, cells[i]->edges[1]->index);
				throw new MazeException(err_string);
			}
		}
		if (cells[i]->edges[2]) {
			if (cells[i]->edges[2]->neighbors[0] == (Cell*)i)
				cells[i]->edges[2]->neighbors[0] = cells[i];
			else if (cells[i]->edges[2]->neighbors[1] == (Cell*)i)
				cells[i]->edges[2]->neighbors[1] = cells[i];
			else {
				sprintf(err_string,
					"Maze: Cell %d not one of edge %d's neighbors",
					i, cells[i]->edges[2]->index);
				throw new MazeException(err_string);
			}
		}
		if (cells[i]->edges[3]) {
			if (cells[i]->edges[3]->neighbors[0] == (Cell*)i)
				cells[i]->edges[3]->neighbors[0] = cells[i];
			else if (cells[i]->edges[3]->neighbors[1] == (Cell*)i)
				cells[i]->edges[3]->neighbors[1] = cells[i];
			else {
				sprintf(err_string,
					"Maze: Cell %d not one of edge %d's neighbors",
					i, cells[i]->edges[3]->index);
				throw new MazeException(err_string);
			}
		}
	}

	if (fscanf(f, "%g %g %g %g %g",
		&(viewer_posn[X]), &(viewer_posn[Y]), &(viewer_posn[Z]),
		&(viewer_dir), &(viewer_fov)) != 5)
		throw new MazeException("Maze: Error reading view information.");

	// Some edges have no neighbor on one side, so be sure to set their
	// pointers to NULL. (They were set at -1 by the save/load process.)
	for (i = 0; i < num_edges; i++) {
		if (edges[i]->neighbors[0] == (Cell*)-1)
			edges[i]->neighbors[0] = NULL;
		if (edges[i]->neighbors[1] == (Cell*)-1)
			edges[i]->neighbors[1] = NULL;
	}

	fclose(f);

	Set_Extents();

	// Figure out which cell the viewer is in, starting off by guessing the
	// 0th cell.
	Find_View_Cell(cells[0]);

	frame_num = 0;
}


//**********************************************************************
//
// * Destructor must free all the memory allocated.
//======================================================================
Maze::
~Maze(void)
//======================================================================
{
	int i;

	for (i = 0; i < num_vertices; i++)
		delete vertices[i];
	delete[] vertices;

	for (i = 0; i < num_edges; i++)
		delete edges[i];
	delete[] edges;

	for (i = 0; i < num_cells; i++)
		delete cells[i];
	delete[] cells;
}


//**********************************************************************
//
// * Randomly generate the edge's opaque and transparency for an empty maze
//======================================================================
void Maze::
Build_Connectivity(const int num_x, const int num_y,
	const float sx, const float sy)
	//======================================================================
{
	int	i, j, k;
	int edge_i;

	// Ugly code to allocate all the memory for a new maze and to associate
	// edges with vertices and faces with edges.

	// Allocate and position the vertices.
	num_vertices = (num_x + 1) * (num_y + 1);
	vertices = new Vertex * [num_vertices];
	k = 0;
	for (i = 0; i < num_y + 1; i++) {
		for (j = 0; j < num_x + 1; j++) {
			vertices[k] = new Vertex(k, j * sx, i * sy);
			k++;
		}
	}

	// Allocate the edges, and associate them with their vertices.
	// Edges in the x direction get the first num_x * ( num_y + 1 ) indices,
	// edges in the y direction get the rest.
	num_edges = (num_x + 1) * num_y + (num_y + 1) * num_x;
	edges = new Edge * [num_edges];
	k = 0;
	for (i = 0; i < num_y + 1; i++) {
		int row = i * (num_x + 1);
		for (j = 0; j < num_x; j++) {
			int vs = row + j;
			int ve = row + j + 1;
			edges[k] = new Edge(k, vertices[vs], vertices[ve],
				rand() / (float)RAND_MAX * 0.5f + 0.25f,
				rand() / (float)RAND_MAX * 0.5f + 0.25f,
				rand() / (float)RAND_MAX * 0.5f + 0.25f);
			k++;
		}
	}

	edge_i = k;
	for (i = 0; i < num_y; i++) {
		int row = i * (num_x + 1);
		for (j = 0; j < num_x + 1; j++) {
			int vs = row + j;
			int ve = row + j + num_x + 1;
			edges[k] = new Edge(k, vertices[vs], vertices[ve],
				rand() / (float)RAND_MAX * 0.5f + 0.25f,
				rand() / (float)RAND_MAX * 0.5f + 0.25f,
				rand() / (float)RAND_MAX * 0.5f + 0.25f);
			k++;
		}
	}

	// Allocate the cells and associate them with their edges.
	num_cells = num_x * num_y;
	cells = new Cell * [num_cells];
	k = 0;
	for (i = 0; i < num_y; i++) {
		int row_x = i * (num_x + 1);
		int row_y = i * num_x;
		for (j = 0; j < num_x; j++) {
			int px = edge_i + row_x + 1 + j;
			int py = row_y + j + num_x;
			int mx = edge_i + row_x + j;
			int my = row_y + j;
			cells[k] = new Cell(k, edges[px], edges[py], edges[mx], edges[my]);
			edges[px]->Add_Cell(cells[k], Edge::LEFT);
			edges[py]->Add_Cell(cells[k], Edge::RIGHT);
			edges[mx]->Add_Cell(cells[k], Edge::RIGHT);
			edges[my]->Add_Cell(cells[k], Edge::LEFT);
			k++;
		}
	}
}


//**********************************************************************
//
// * Add edges from cell to the set that are available for removal to
//   grow the maze.
//======================================================================
static void
Add_To_Available(Cell* cell, int* available, int& num_available)
//======================================================================
{
	int i, j;

	// Add edges from cell to the set that are available for removal to
	// grow the maze.

	for (i = 0; i < 4; i++) {
		Cell* neighbor = cell->edges[i]->Neighbor(cell);

		if (neighbor && !neighbor->counter) {
			int candidate = cell->edges[i]->index;
			for (j = 0; j < num_available; j++)
				if (candidate == available[j]) {
					printf("Breaking early\n");
					break;
				}
			if (j == num_available) {
				available[num_available] = candidate;
				num_available++;
			}
		}
	}

	cell->counter = 1;
}


//**********************************************************************
//
// * Grow a maze by removing candidate edges until all the cells are
//   connected. The edges are not actually removed, they are just made
//   transparent.
//======================================================================
void Maze::
Build_Maze()
//======================================================================
{
	Cell* to_expand;
	int     index;
	int* available = new int[num_edges];
	int     num_available = 0;
	int	    num_visited;
	int	    i;

	srand(time(NULL));

	// Choose a random starting cell.
	index = (int)floor((rand() / (float)RAND_MAX) * num_cells);
	to_expand = cells[index];
	Add_To_Available(to_expand, available, num_available);
	num_visited = 1;

	// Join cells up by making edges opaque.
	while (num_visited < num_cells && num_available > 0) {
		int ei;

		index = (int)floor((rand() / (float)RAND_MAX) * num_available);
		to_expand = NULL;

		ei = available[index];

		if (edges[ei]->neighbors[0] &&
			!edges[ei]->neighbors[0]->counter)
			to_expand = edges[ei]->neighbors[0];
		else if (edges[ei]->neighbors[1] &&
			!edges[ei]->neighbors[1]->counter)
			to_expand = edges[ei]->neighbors[1];

		if (to_expand) {
			edges[ei]->opaque = false;
			Add_To_Available(to_expand, available, num_available);
			num_visited++;
		}

		available[index] = available[num_available - 1];
		num_available--;
	}

	for (i = 0; i < num_cells; i++)
		cells[i]->counter = 0;
}


//**********************************************************************
//
// * Go through all the vertices looking for the minimum and maximum
//   extents of the maze.
//======================================================================
void Maze::
Set_Extents(void)
//======================================================================
{
	int i;

	min_xp = vertices[0]->posn[Vertex::X];
	max_xp = vertices[0]->posn[Vertex::X];
	min_yp = vertices[0]->posn[Vertex::Y];
	max_yp = vertices[0]->posn[Vertex::Y];
	for (i = 1; i < num_vertices; i++) {
		if (vertices[i]->posn[Vertex::X] > max_xp)
			max_xp = vertices[i]->posn[Vertex::X];
		if (vertices[i]->posn[Vertex::X] < min_xp)
			min_xp = vertices[i]->posn[Vertex::X];
		if (vertices[i]->posn[Vertex::Y] > max_yp)
			max_yp = vertices[i]->posn[Vertex::Y];
		if (vertices[i]->posn[Vertex::Y] < min_yp)
			min_yp = vertices[i]->posn[Vertex::Y];
	}
}


//**********************************************************************
//
// * Figure out which cell the view is in, using seed_cell as an
//   initial guess. This procedure works by repeatedly checking
//   whether the viewpoint is in the current cell. If it is, we're
//   done. If not, Point_In_Cell returns in new_cell the next cell
//   to test. The new cell is the one on the other side of an edge
//   that the point is "outside" (meaning that it might be inside the
//   new cell).
//======================================================================
void Maze::
Find_View_Cell(Cell* seed_cell)
//======================================================================
{
	Cell* new_cell;

	// 
	while (!(seed_cell->Point_In_Cell(viewer_posn[X], viewer_posn[Y],
		viewer_posn[Z], new_cell))) {
		if (new_cell == 0) {
			// The viewer is outside the top or bottom of the maze.
			throw new MazeException("Maze: View not in maze\n");
		}

		seed_cell = new_cell;
	}

	view_cell = seed_cell;
}


//**********************************************************************
//
// * Move the viewer's position. This method will do collision detection
//   between the viewer's location and the walls of the maze and prevent
//   the viewer from passing through walls.
//======================================================================
void Maze::
Move_View_Posn(const float dx, const float dy, const float dz)
//======================================================================
{
	Cell* new_cell;
	float   xs, ys, zs, xe, ye, ze;

	// Move the viewer by the given amount. This does collision testing to
	// prevent walking through walls. It also keeps track of which cells the
	// viewer is in.

	// Set up a line segment from the start to end points of the motion.
	xs = viewer_posn[X];
	ys = viewer_posn[Y];
	zs = viewer_posn[Z];
	xe = xs + dx;
	ye = ys + dy;
	ze = zs + dz;

	// Fix the z to keep it in the maze.
	if (ze > 1.0f - BUFFER)
		ze = 1.0f - BUFFER;
	if (ze < BUFFER - 1.0f)
		ze = BUFFER - 1.0f;

	// Clip_To_Cell clips the motion segment to the view_cell if the
	// segment intersects an opaque edge. If the segment intersects
	// a transparent edge (through which it can pass), then it clips
	// the motion segment so that it _starts_ at the transparent edge,
	// and it returns the cell the viewer is entering. We keep going
	// until Clip_To_Cell returns NULL, meaning we've done as much of
	// the motion as is possible without passing through walls.
	while ((new_cell = view_cell->Clip_To_Cell(xs, ys, xe, ye, BUFFER)))
		view_cell = new_cell;

	// The viewer is at the end of the motion segment, which may have
	// been clipped.
	viewer_posn[X] = xe;
	viewer_posn[Y] = ye;
	viewer_posn[Z] = ze;
}

//**********************************************************************
//
// * Set the viewer's location 
//======================================================================
void Maze::
Set_View_Posn(float x, float y, float z)
//======================================================================
{
	// First make sure it's in some cell.
	// This assumes that the maze is rectangular.
	if (x < min_xp + BUFFER)
		x = min_xp + BUFFER;
	if (x > max_xp - BUFFER)
		x = max_xp - BUFFER;
	if (y < min_yp + BUFFER)
		y = min_yp + BUFFER;
	if (y > max_yp - BUFFER)
		y = max_yp - BUFFER;
	if (z < -1.0f + BUFFER)
		z = -1.0f + BUFFER;
	if (z > 1.0f - BUFFER)
		z = 1.0f - BUFFER;

	viewer_posn[X] = x;
	viewer_posn[Y] = y;
	viewer_posn[Z] = z;

	// Figure out which cell we're in.
	Find_View_Cell(cells[0]);
}


//**********************************************************************
//
// * Set the angle in which the viewer is looking.
//======================================================================
void Maze::
Set_View_Dir(const float d)
//======================================================================
{
	viewer_dir = d;
}


//**********************************************************************
//
// * Set the horizontal field of view.
//======================================================================
void Maze::
Set_View_FOV(const float f)
//======================================================================
{
	viewer_fov = f;
}


//**********************************************************************
//
// * Draws the map view of the maze. It is passed the minimum and maximum
//   corners of the window in which to draw.
//======================================================================
void Maze::
Draw_Map(int min_x, int min_y, int max_x, int max_y)
//======================================================================
{
	int	    height;
	float   scale_x, scale_y, scale;
	int	    i;

	// Figure out scaling factors and the effective height of the window.
	scale_x = (max_x - min_x - 10) / (max_xp - min_xp);
	scale_y = (max_y - min_y - 10) / (max_yp - min_yp);
	scale = scale_x > scale_y ? scale_y : scale_x;
	height = (int)ceil(scale * (max_yp - min_yp));

	min_x += 5;
	min_y += 5;

	// Draw all the opaque edges.
	for (i = 0; i < num_edges; i++)
		if (edges[i]->opaque) {
			float   x1, y1, x2, y2;

			x1 = edges[i]->endpoints[Edge::START]->posn[Vertex::X];
			y1 = edges[i]->endpoints[Edge::START]->posn[Vertex::Y];
			x2 = edges[i]->endpoints[Edge::END]->posn[Vertex::X];
			y2 = edges[i]->endpoints[Edge::END]->posn[Vertex::Y];

			fl_color((unsigned char)floor(edges[i]->color[0] * 255.0),
				(unsigned char)floor(edges[i]->color[1] * 255.0),
				(unsigned char)floor(edges[i]->color[2] * 255.0));
			fl_line_style(FL_SOLID);
			fl_line(min_x + (int)floor((x1 - min_xp) * scale),
				min_y + height - (int)floor((y1 - min_yp) * scale),
				min_x + (int)floor((x2 - min_xp) * scale),
				min_y + height - (int)floor((y2 - min_yp) * scale));
		}
	/*else
	{
		float   x1, y1, x2, y2;

		x1 = edges[i]->endpoints[Edge::START]->posn[Vertex::X];
		y1 = edges[i]->endpoints[Edge::START]->posn[Vertex::Y];
		x2 = edges[i]->endpoints[Edge::END]->posn[Vertex::X];
		y2 = edges[i]->endpoints[Edge::END]->posn[Vertex::Y];

		fl_color((unsigned char)0,
			(unsigned char)0,
			(unsigned char)0);
		fl_line_style(FL_SOLID);
		fl_line(min_x + (int)floor((x1 - min_xp) * scale),
			min_y + height - (int)floor((y1 - min_yp) * scale),
			min_x + (int)floor((x2 - min_xp) * scale),
			min_y + height - (int)floor((y2 - min_yp) * scale));
	}*/
}


//**********************************************************************
//
// * Draws the first-person view of the maze. It is passed the focal distance.
//   THIS IS THE FUINCTION YOU SHOULD MODIFY.
//======================================================================
void Maze::
Draw_View(const float focal_dist, const double m[])
//======================================================================
{
	for (int i = 0; i < 16; i++)
	{
		multiMatrix[i] = m[i];
	}
	frame_num++;
	HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);
	//###################################################################
	// TODO
	// The rest is up to you!
	//###################################################################

	//GL method
	/*glClear(GL_DEPTH_BUFFER_BIT);

	glEnable(GL_DEPTH_TEST);*/

	//printf("view: ( %f, %f, %f),dir: %f\n", viewer_posn[0], viewer_posn[1], viewer_posn[2], viewer_dir);
	//for (int i = 0; i < (int)this->num_edges; i++)
	//{
	//	float color[3] = { this->edges[i]->color[0],this->edges[i]->color[1], this->edges[i]->color[2] };//rgb
	//	float edge_start[2] = { this->edges[i]->endpoints[Edge::START]->posn[Vertex::X], this->edges[i]->endpoints[Edge::START]->posn[Vertex::Y] };
	//	float edge_end[2] = { this->edges[i]->endpoints[Edge::END]->posn[Vertex::X], this->edges[i]->endpoints[Edge::END]->posn[Vertex::Y] };

	//	//printf("( %f, %f)->( %f, %f),(%f,%f,%f)\n",edge_start[0], edge_start[1], edge_end[0], edge_end[1], 255*color[0], 255*color[1], 255*color[2]);
	//	if (this->edges[i]->opaque)//can see
	//	{
	//		if (i == 3)
	//		{
	//			p = true;
	//		}
	//		else
	//		{
	//			p = false;
	//		}
	//		Draw_Wall(edge_start, edge_end, color);
	//	}
	//}
	for (int i = 0; i < num_edges; i++)
	{
		isUsed[i] = -1;
		isTransperant[i] = -1;
	}
	/*int index = 0;
	for (int i = 0; i < num_edges; i++)
	{
		if (!this->edges[i]->opaque)
		{
			isUsed[index] = i;
			index++;
		}
	}*/
	isTransperantIndex = 0;
	max_radius = 1;
	double cameraPos[2] = { viewer_posn[0], viewer_posn[1] };
	//(double) middle[2] = { max_radius*cos(viewer_dir * RAG_Div), max_radius*sin(viewer_dir * RAG_Div) };
	double leftClip[2] = { max_radius * cos(To_Radians(viewer_dir + viewer_fov / 2)),max_radius * sin(To_Radians(viewer_dir + viewer_fov / 2)) };

	double rightClip[2] = { max_radius * cos(To_Radians(viewer_dir - viewer_fov / 2)),max_radius * sin(To_Radians(viewer_dir - viewer_fov / 2)) };
	if (Display)
	{
		system("cls");
		SetColor();
		printf("DRAW START!\n");
	}
	Draw_Cell(cameraPos, leftClip, rightClip);
}

//**********************************************************************
//
// * check which cell we need to draw
//   THIS IS THE FUINCTION YOU SHOULD MODIFY.
//======================================================================
void Maze::
Draw_Cell(const double pos[], double left[], double right[])
//======================================================================
{
	//printf("Pos: ( %f, %f), L: ( %f, %f), R: ( %f, %f)\n", pos[0], pos[1], left[0], left[1], right[0], right[1]);

	int count = 0;
	drawCoordinate endCoor = Find_Edge(pos[0], pos[1], pos[0] + right[0], pos[1] + right[1], false);//Find Right Wall
	double RClip[2] = { pos[0] + right[0] * endCoor.crossIn, pos[1] + right[1] * endCoor.crossIn };//the Right screen coordinate
	right[0] = RClip[0] - pos[0];//the right vertex
	right[1] = RClip[1] - pos[1];//the right vertex

	drawCoordinate startCoor = Find_Edge(pos[0], pos[1], pos[0] + left[0], pos[1] + left[1], false);//Find Left Wall
	double LClip[2] = { pos[0] + left[0] * startCoor.crossIn,pos[1] + left[1] * startCoor.crossIn };//the Light screen coordinate
	if (Display)
	{
		SetColor(11);
		printf("Pos: ( %f, %f), L: ( %f, %f), R: ( %f, %f)\n", pos[0], pos[1], left[0], left[1], right[0], right[1]);
		printf("L:%d: (%f, %f)\n", startCoor.index, LClip[0], LClip[1]);
		printf("R:%d: (%f, %f)\n", endCoor.index, RClip[0], RClip[1]);
		SetColor();
	}
	double cosine = (left[0] * right[0]) + (left[1] * right[1]);
	cosine /= sqrt(pow(left[0], 2) + pow(left[1], 2)); cosine /= sqrt(pow(right[0], 2) + pow(right[1], 2));
	while ((fabs(cosine - 1) > 0.0000001))
	{
		double DrawClipPoint[2][2] = {
		RClip[0],RClip[1],
		LClip[0],LClip[1],
		};
		count++;
		if (count > 10)
		{
			if (Display)
			{
				SetColor(15);
				printf("wrong\n");
				SetColor();
			}

			return;
		}
		drawCoordinate edgeIndex = Find_Edge(pos[0], pos[1], LClip[0], LClip[1], true);//Find Left Wall
		double tempLeft[2] = { left[0] * cos(To_Radians(-0.000001))-left[1]* sin(To_Radians(-0.000001)), left[0] * sin(To_Radians(-0.000001)) + left[1] * cos(To_Radians(-0.000001)) };
		drawCoordinate tempEdgeIndex = Find_Edge(pos[0], pos[1], pos[0] + tempLeft[0], pos[1] + tempLeft[1], false);//Find Left Wall
		if ((tempEdgeIndex.index != -1) && (tempEdgeIndex.crossIn < edgeIndex.crossIn))
		{
			edgeIndex = tempEdgeIndex;
		}
		if (edgeIndex.index != -1)
		{
			if (this->edges[edgeIndex.index]->opaque)
			{
				if (Display)
				{
					SetColor(14);
					printf("index: %d, opaque: ", edgeIndex.index);
				}

				double edge_draw[2][2] =
				{
					{ my_edge[edgeIndex.index][0], my_edge[edgeIndex.index][1] },//start
					{ my_edge[edgeIndex.index][2], my_edge[edgeIndex.index][3] }
				};
				//printf("1: start: ( %f, %f ), end: ( %f , %f )\n", edge_draw[0][0], edge_draw[0][1], edge_draw[1][0], edge_draw[1][1]);
				float color[3] = { (float)my_edge[edgeIndex.index][4], (float)my_edge[edgeIndex.index][5], (float)my_edge[edgeIndex.index][6] };//rgb
				Draw_Part(edge_draw, pos, left, right, false);
				float edge_start[2] = { edge_draw[0][0], edge_draw[0][1] };
				float edge_end[2] = { edge_draw[1][0], edge_draw[1][1] };
				if (Display)
				{
					printf("4: start: ( %f, %f ), end: ( %f , %f )\n", edge_draw[0][0], edge_draw[0][1], edge_draw[1][0], edge_draw[1][1]);
				}
				Draw_Wall(edge_start, edge_end, color);
				if (Display)
				{
					Draw_Point(DrawClipPoint[0], false);
					Draw_Point(DrawClipPoint[1], true);
					Draw_Point(edge_draw[0], false);
					Draw_Point(edge_draw[1], true);
				}
				
				if (endCoor.index == edgeIndex.index)
				{
					LClip[0] = RClip[0];
					LClip[1] = RClip[1];
					left[0] = LClip[0] - pos[0];
					left[1] = LClip[1] - pos[1];
				}
				else
				{
					LClip[0] = edge_draw[1][0];
					LClip[1] = edge_draw[1][1];
					left[0] = (1000 * LClip[0] - 1000 * (double)pos[0]) / 1000;
					left[1] = (1000 * LClip[1] - 1000 * (double)pos[1]) / 1000;
					if (Display)
					{
						printf("( %f, %f)\n", LClip[0], LClip[1]);
					}
				}
			}
			else
			{
				//system("pause");
				if (Display)
				{
					SetColor(10);
					printf("index: %d, transparent: ", edgeIndex.index);
				}
				double edge_throw[2][2] =
				{
					{ my_edge[edgeIndex.index][0], my_edge[edgeIndex.index][1] },//start
					{ my_edge[edgeIndex.index][2], my_edge[edgeIndex.index][3] }
				};

				isTransperant[isTransperantIndex] = edgeIndex.index;
				isTransperantIndex++;


				Draw_Part(edge_throw, pos, left, right, true);
				//Draw_Part(edge_throw, pos, left, right, false);
				float edge_start[2] = { edge_throw[0][0], edge_throw[0][1] };
				float edge_end[2] = { edge_throw[1][0], edge_throw[1][1] };
				if (Display)
				{
					Draw_Point(DrawClipPoint[0], false);
					Draw_Point(DrawClipPoint[1], true);
					Draw_Point(edge_throw[0], false);
					Draw_Point(edge_throw[1], true);
				}
				left[0] = (1000 * edge_throw[0][0] - 1000 * (double)pos[0]) / 1000;
				left[1] = (1000 * edge_throw[0][1] - 1000 * (double)pos[1]) / 1000;


				double rightTemp[2] = { (1000 * edge_throw[1][0] - 1000 * (double)pos[0]) / 1000 , (1000 * edge_throw[1][1] - 1000 * (double)pos[1]) / 1000 };

				//drawCoordinate RedgeIndex;
				//if (find)
				//{
				//	isTransperantIndex++;
				//	RedgeIndex = Find_Edge(pos[0], pos[1], edge_throw[1][0], edge_throw[1][1], false);//Find Left Wall
				//	if (RedgeIndex.index != -1)
				//	{
				//		rightTemp[0] = RedgeIndex.crossIn * rightTemp[0];
				//		rightTemp[1] = RedgeIndex.crossIn * rightTemp[1];
				//		printf("index: %d, ", RedgeIndex.index);
				//	}
				//	isTransperantIndex--;
				//}
				//printf("left: (%f, %f), rightTemp: (%f, %f)\n", left[0], left[1], rightTemp[0], rightTemp[1]);

				if (Display)
				{
					printf("l: (%f, %f), r: (%f, %f)\n", pos[0] + left[0], pos[1] + left[1], pos[0] + rightTemp[0], pos[1] + rightTemp[1]);
				}
				isTransperantIndex++;
				Draw_Cell(pos, left, rightTemp);
				LClip[0] = edge_throw[1][0];
				LClip[1] = edge_throw[1][1];
				left[0] = rightTemp[0];
				left[1] = rightTemp[1];
			}
			//printf("left: (%f, %f), right: (%f, %f)\n", left[0], left[1], right[0], right[1]);
			if (Display)
			{
				SetColor();
				printf("LClip: (%f, %f), RClip: (%f, %f)\n", LClip[0], LClip[1], RClip[0], RClip[1]);
			}
		}
		cosine = (left[0] * right[0]) + (left[1] * right[1]);
		cosine /= sqrt(pow(left[0], 2) + pow(left[1], 2)); cosine /= sqrt(pow(right[0], 2) + pow(right[1], 2));
		if (Display)
		{
			SetColor();
			printf("cosine: %f\n", cosine);
		}
	}
	if (Display)
	{
		SetColor(15);
		printf("success\n");
		SetColor();
	}
	return;
}

//**********************************************************************
//
// * find the nearest edge, and return its index
//   THIS IS THE FUINCTION YOU SHOULD MODIFY.
//======================================================================
drawCoordinate	Maze::
Find_Edge(double xs, double ys, double xe, double ye, bool Clip)
{
	drawCoordinate coordinate;
	coordinate.index = -1;
	double distance = 100000;
	LineSeg in_seg(xs, ys, xe, ye);
	for (int i = 0; i < (int)this->num_edges; i++)
	{
		bool k = false;
		bool l = false;
		bool limit = false;
		//if (Clip)
		//{
			for (int j = 0; j < (int)this->num_edges; j++)
			{
				if (isUsed[j] == i)
				{
					k = true;
					//printf("true\n");
				}
			}
		//}
		/*else if (isTransperant[isTransperantIndex - 1] == -1)
		{
			if (((fabs(xe - this->my_edge[i][0]) > 0.5) || (fabs(ye - this->my_edge[i][1]) > 0.5)) && ((fabs(xe - this->my_edge[i][2]) > 0.5) || (fabs(ye - this->my_edge[i][3]) > 0.5)))
			{
				limit = false;

			}
			else
			{
				limit = true;
			}
		}*/
		for (int j = 0; j < isTransperantIndex - 1; j++)
		{
			if (isTransperant[j] == i)
			{
				l = true;
				//printf("true\n");
			}
		}


		if (!k && !l && !limit)
		{
			if (Clip)
			{
				//printf("false: %d\n", i);
			}
			LineSeg e(this->my_edge[i][0], this->my_edge[i][1], this->my_edge[i][2], this->my_edge[i][3]);
			double crossIn = in_seg.Cross_Param(e);
			double crossE = e.Cross_Param(in_seg);
			if (Clip)
			{
				//printf("cross: %f, %f\n", crossIn, crossE);
			}
			if ((crossIn >= 0) && (crossE >= 0) && crossE <= 1.0001)
			{
				if (crossIn < distance)
				{
					coordinate.index = i;
					coordinate.crossIn = crossIn;
					coordinate.crossE = crossE;
					distance = crossIn;
					if (Clip)
					{
						//printf("min\n");
						//printf("%d ", i);
					}
				}
			}
		}
	}

	if (Clip)
	{
		int i = 0;
		while (isUsed[i] != -1)
		{
			i++;
		}
		isUsed[i] = coordinate.index;
		//printf("\n");
	}
	return coordinate;
}


// check which part of wall need to draw
// return edge[0]:left, edge[1]:right
void	Maze
::Draw_Part(double edge[][2], const double pos[], double left[], double right[], bool isTrans)
{
	// left X edgeTemp -> crossVector 
	// if((0,0,-1)) don't change, else start = end end = start
	// if(a1b2-a2b1)<0
	double edgeTemp[2] = { edge[1][0] - edge[0][0], edge[1][1] - edge[0][1] };//end-start
	double crossVector = left[0] * edgeTemp[1] - left[1] * edgeTemp[0];
	if (crossVector > 0)
	{
		double temp[2] = { edge[0][0],edge[0][1] };
		edge[0][0] = edge[1][0];
		edge[0][1] = edge[1][1];
		edge[1][0] = temp[0];
		edge[1][1] = temp[1];
	}
	//printf("2: start: ( %f, %f ), end: ( %f, %f )\n", edge[0][0], edge[0][1], edge[1][0], edge[1][1]);


	double start[2] = { edge[0][0], edge[0][1] };
	edgeTemp[0] = edge[1][0] - edge[0][0]; edgeTemp[1] = edge[1][1] - edge[0][1];

	LineSeg wall_seg(edge[0][0], edge[0][1], edge[1][0], edge[1][1]);
	LineSeg left_seg(pos[0], pos[1], pos[0] + left[0], pos[1] + left[1]);
	LineSeg right_seg(pos[0], pos[1], pos[0] + right[0], pos[1] + right[1]);
	double crossL = wall_seg.Cross_Param(left_seg);
	double crossR = wall_seg.Cross_Param(right_seg);
	double RcrossWall = right_seg.Cross_Param(wall_seg);
	if (isTrans)
	{
		crossR -= 0.0001;
	}
	edge[0][0] = start[0] + crossL * edgeTemp[0];
	edge[0][1] = start[1] + crossL * edgeTemp[1];
	if ((RcrossWall >= 0) && (crossR >= 0) && (crossR <= 1.00001))
	{
		edge[1][0] = start[0] + crossR * edgeTemp[0];
		edge[1][1] = start[1] + crossR * edgeTemp[1];
	}

	//printf("3: crossR: %f start: ( %f, %f ), end: ( %f, %f )\n", crossR, edge[0][0], edge[0][1], edge[1][0], edge[1][1]);
}


void Maze
::Draw_Point(const double dot[2], const bool L) {
	glMatrixMode(GL_PROJECTION);//camera
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	double startTemp0[4] = { dot[Y], 0,dot[X] ,1 };
	double edgeStart0[4] = { 0 };
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if ((multiMatrix[i * 4 + j] != 0.0f) && (startTemp0[j] != 0.0f))
				edgeStart0[i] += multiMatrix[i * 4 + j] * startTemp0[j];
		}
	}

	double temp = -0.05;
	if (L)
	{
		temp += 0.05;
	}

	glPointSize(10);
	glBegin(GL_POINTS);
	glColor3f(1.0f, 1.0f, 1.0f);
	glVertex4f(edgeStart0[0] + temp, edgeStart0[1], edgeStart0[2], edgeStart0[3]);
	glEnd();
}

//**********************************************************************
//
// *  Draws the wall of the maze.
//   THIS IS THE FUINCTION YOU SHOULD MODIFY.
//======================================================================
void Maze::
Draw_Wall(const float start[2], const float end[2], const float color[3])
//======================================================================
{
	glMatrixMode(GL_PROJECTION);//camera
	glLoadIdentity();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	double edgeTemp[4][4] = {
		{ start[Y], 1.0f,start[X] ,1 },//start0
		{ end[Y], 1.0f,end[X] ,1 },//end0
		{ end[Y], -1.0f,end[X] ,1 },//end1
		{ start[Y], -1.0f,start[X] ,1 }//start1
	};
	double edgeDraw[4][4] = {
		{ 0,0,0,0 },//start0
		{ 0,0,0,0 },//end0
		{ 0,0,0,0 },//end1
		{ 0,0,0,0 }//start1
	};
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				if ((multiMatrix[j * 4 + k] != 0.0f) && (edgeTemp[i][k] != 0.0f))
					edgeDraw[i][j] += multiMatrix[j * 4 + k] * edgeTemp[i][k];
			}
			//if(p)
			//printf("%f ", edgeDraw[i][j]);
		}
		//if (p)
		//printf("\n");
	}

	glBegin(GL_POLYGON);
	glColor3fv(color);
	for (int i = 0; i < 4; i++)
	{
		if (edgeDraw[i][2] > 0)
		{
			glVertex2f(edgeDraw[i][0] / edgeDraw[i][2], edgeDraw[i][1] / edgeDraw[i][2]);
		}
	}
	/*glVertex4f(edgeStart0[0], edgeStart0[1], edgeStart0[2], edgeStart0[3]);
	glVertex4f(edgeEnd0[0], edgeEnd0[1], edgeEnd0[2], edgeEnd0[3]);
	glVertex4f(edgeEnd1[0], edgeEnd1[1], edgeEnd1[2], edgeEnd1[3]);
	glVertex4f(edgeStart1[0], edgeStart1[1], edgeStart1[2], edgeStart1[3]);*/
	glEnd();
}


//**********************************************************************
//
// * Draws the frustum on the map view of the maze. It is passed the
//   minimum and maximum corners of the window in which to draw.
//======================================================================
void Maze::
Draw_Frustum(int min_x, int min_y, int max_x, int max_y)
//======================================================================
{
	int	  height;
	float   scale_x, scale_y, scale;
	float   view_x, view_y;

	// Draws the view frustum in the map. Sets up all the same viewing
	// parameters as draw().
	scale_x = (max_x - min_x - 10) / (max_xp - min_xp);
	scale_y = (max_y - min_y - 10) / (max_yp - min_yp);
	scale = scale_x > scale_y ? scale_y : scale_x;
	height = (int)ceil(scale * (max_yp - min_yp));

	min_x += 5;
	min_y += 5;

	view_x = (viewer_posn[X] - min_xp) * scale;
	view_y = (viewer_posn[Y] - min_yp) * scale;
	fl_line(min_x + (int)floor(view_x +
		cos(To_Radians(viewer_dir + viewer_fov / 2.0)) * scale),
		min_y + height -
		(int)floor(view_y +
			sin(To_Radians(viewer_dir + viewer_fov / 2.0)) *
			scale),
		min_x + (int)floor(view_x),
		min_y + height - (int)floor(view_y));
	fl_line(min_x + (int)floor(view_x +
		cos(To_Radians(viewer_dir - viewer_fov / 2.0)) *
		scale),
		min_y + height -
		(int)floor(view_y + sin(To_Radians(viewer_dir - viewer_fov / 2.0)) *
			scale),
		min_x + (int)floor(view_x),
		min_y + height - (int)floor(view_y));
}


//**********************************************************************
//
// * Draws the viewer's cell and its neighbors in the map view of the maze.
//   It is passed the minimum and maximum corners of the window in which
//   to draw.
//======================================================================
void Maze::
Draw_Neighbors(int min_x, int min_y, int max_x, int max_y)
//======================================================================
{
	int	    height;
	float   scale_x, scale_y, scale;
	int	    i, j;

	// Draws the view cell and its neighbors in the map. This works
	// by drawing just the neighbor's edges if there is a neighbor,
	// otherwise drawing the edge. Every edge is shared, so drawing the
	// neighbors' edges also draws the view cell's edges.

	scale_x = (max_x - min_x - 10) / (max_xp - min_xp);
	scale_y = (max_y - min_y - 10) / (max_yp - min_yp);
	scale = scale_x > scale_y ? scale_y : scale_x;
	height = (int)ceil(scale * (max_yp - min_yp));

	min_x += 5;
	min_y += 5;

	for (i = 0; i < 4; i++) {
		Cell* neighbor = view_cell->edges[i]->Neighbor(view_cell);

		if (neighbor) {
			for (j = 0; j < 4; j++) {
				Edge* e = neighbor->edges[j];

				if (e->opaque) {
					float   x1, y1, x2, y2;

					x1 = e->endpoints[Edge::START]->posn[Vertex::X];
					y1 = e->endpoints[Edge::START]->posn[Vertex::Y];
					x2 = e->endpoints[Edge::END]->posn[Vertex::X];
					y2 = e->endpoints[Edge::END]->posn[Vertex::Y];

					fl_color((unsigned char)floor(e->color[0] * 255.0),
						(unsigned char)floor(e->color[1] * 255.0),
						(unsigned char)floor(e->color[2] * 255.0));
					fl_line_style(FL_SOLID);
					fl_line(min_x + (int)floor((x1 - min_xp) * scale),
						min_y + height - (int)floor((y1 - min_yp) * scale),
						min_x + (int)floor((x2 - min_xp) * scale),
						min_y + height - (int)floor((y2 - min_yp) * scale));
				}
			}
		}
		else {
			Edge* e = view_cell->edges[i];

			if (e->opaque) {
				float   x1, y1, x2, y2;

				x1 = e->endpoints[Edge::START]->posn[Vertex::X];
				y1 = e->endpoints[Edge::START]->posn[Vertex::Y];
				x2 = e->endpoints[Edge::END]->posn[Vertex::X];
				y2 = e->endpoints[Edge::END]->posn[Vertex::Y];

				fl_color((unsigned char)floor(e->color[0] * 255.0),
					(unsigned char)floor(e->color[1] * 255.0),
					(unsigned char)floor(e->color[2] * 255.0));
				fl_line_style(FL_SOLID);
				fl_line(min_x + (int)floor((x1 - min_xp) * scale),
					min_y + height - (int)floor((y1 - min_yp) * scale),
					min_x + (int)floor((x2 - min_xp) * scale),
					min_y + height - (int)floor((y2 - min_yp) * scale));
			}
		}
	}
}


//**********************************************************************
//
// * Save the maze to a file of the given name.
//======================================================================
bool Maze::
Save(const char* filename)
//======================================================================
{
	FILE* f = fopen(filename, "w");
	int	    i;

	// Dump everything to a file of the given name. Returns false if it
	// couldn't open the file. True otherwise.

	if (!f) {
		return false;
	}

	fprintf(f, "%d\n", num_vertices);
	for (i = 0; i < num_vertices; i++)
		fprintf(f, "%g %g\n", vertices[i]->posn[Vertex::X],
			vertices[i]->posn[Vertex::Y]);

	fprintf(f, "%d\n", num_edges);
	for (i = 0; i < num_edges; i++)
		fprintf(f, "%d %d %d %d %d %g %g %g\n",
			edges[i]->endpoints[Edge::START]->index,
			edges[i]->endpoints[Edge::END]->index,
			edges[i]->neighbors[Edge::LEFT] ?
			edges[i]->neighbors[Edge::LEFT]->index : -1,
			edges[i]->neighbors[Edge::RIGHT] ?
			edges[i]->neighbors[Edge::RIGHT]->index : -1,
			edges[i]->opaque ? 1 : 0,
			edges[i]->color[0], edges[i]->color[1], edges[i]->color[2]);

	fprintf(f, "%d\n", num_cells);
	for (i = 0; i < num_cells; i++)
		fprintf(f, "%d %d %d %d\n",
			cells[i]->edges[0] ? cells[i]->edges[0]->index : -1,
			cells[i]->edges[1] ? cells[i]->edges[1]->index : -1,
			cells[i]->edges[2] ? cells[i]->edges[2]->index : -1,
			cells[i]->edges[3] ? cells[i]->edges[3]->index : -1);

	fprintf(f, "%g %g %g %g %g\n",
		viewer_posn[X], viewer_posn[Y], viewer_posn[Z],
		viewer_dir, viewer_fov);

	fclose(f);

	return true;
}