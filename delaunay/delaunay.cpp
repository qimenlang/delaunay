#include "delaunay.h"

namespace dt {

template<typename T>
const std::vector<typename Delaunay<T>::TriangleType>&
Delaunay<T>::triangulate(std::vector<VertexType> &vertices)
{
	// Store the vertices locally
	_vertices = vertices;

	// Determinate the super triangle
	T minX = vertices[0].x;
	T minY = vertices[0].y;
	T maxX = minX;
	T maxY = minY;

	for(std::size_t i = 0; i < vertices.size(); ++i)
	{
		if (vertices[i].x < minX) minX = vertices[i].x;
		if (vertices[i].y < minY) minY = vertices[i].y;
		if (vertices[i].x > maxX) maxX = vertices[i].x;
		if (vertices[i].y > maxY) maxY = vertices[i].y;
	}

	const T dx = maxX - minX;
	const T dy = maxY - minY;
	const T deltaMax = std::max(dx, dy);
	const T midx = (minX + maxX) / 2;
	const T midy = (minY + maxY) / 2;

	const VertexType p1(midx - 20 * deltaMax, midy - deltaMax);
	const VertexType p2(midx, midy + 20 * deltaMax);
	const VertexType p3(midx + 20 * deltaMax, midy - deltaMax);

	// Create a list of triangles, and add the supertriangle in it
	_triangles.push_back(TriangleType(p1, p2, p3));

	for(auto p = begin(vertices); p != end(vertices); p++)
	{
		std::vector<EdgeType> polygon;

		for(auto & t : _triangles)
		{
			if(t.circumCircleContains(*p))
			{
				t.isBad = true;
				polygon.push_back(Edge<T>{*t.a, *t.b});
				polygon.push_back(Edge<T>{*t.b, *t.c});
				polygon.push_back(Edge<T>{*t.c, *t.a});
			}
		}

		_triangles.erase(std::remove_if(begin(_triangles), end(_triangles), [](TriangleType &t){
			return t.isBad;
		}), end(_triangles));

		for(auto e1 = begin(polygon); e1 != end(polygon); ++e1)
		{
			for(auto e2 = e1 + 1; e2 != end(polygon); ++e2)
			{
				if(almost_equal(*e1, *e2))
				{
					e1->isBad = true;
					e2->isBad = true;
				}
			}
		}

		polygon.erase(std::remove_if(begin(polygon), end(polygon), [](EdgeType &e){
			return e.isBad;
		}), end(polygon));

		for(const auto e : polygon)
			_triangles.push_back(TriangleType(*e.v, *e.w, *p));

	}

	_triangles.erase(std::remove_if(begin(_triangles), end(_triangles), [p1, p2, p3](TriangleType &t){
		return t.containsVertex(p1) || t.containsVertex(p2) || t.containsVertex(p3);
	}), end(_triangles));

	for(const auto t : _triangles)
	{
		_edges.push_back(Edge<T>{*t.a, *t.b});
		_edges.push_back(Edge<T>{*t.b, *t.c});
		_edges.push_back(Edge<T>{*t.c, *t.a});
	}

	return _triangles;
}

template<typename T>
std::vector<VertexWeight> Delaunay<T>::getBrayCoord(const VertexType& sample_point) {

	std::vector<VertexWeight> result_weights = { {-1,0.0},{-1,0.0},{-1,0.0} };
	//遍历所有三角形，如果采样点在三角形内部，计算并返回重心坐标
	std::array<T, 3> bary_coord;
	std::array<int, 3> vertex_index;
	for (int i = 0; i < _triangles.size(); i++)
	{
		bary_coord = _triangles[i].GetBaryCentricCoord(sample_point);
		if (bary_coord[0] > 0 && bary_coord[1] > 0 && bary_coord[2] > 0)
		{
			std::cout <<"Triangle Bray Index:"<< i << ":(" << bary_coord[0] << "," << bary_coord[1] << "," << bary_coord[2] << ")" << std::endl;
			for (int j = 0; j < _vertices.size(); j++) {
				if (almost_equal(*_triangles[i].a, _vertices[j])) result_weights[0].vertex_index = j;
				if (almost_equal(*_triangles[i].b, _vertices[j])) result_weights[1].vertex_index = j;
				if (almost_equal(*_triangles[i].c, _vertices[j])) result_weights[2].vertex_index = j;
			}

			for (int k = 0; k < 3; k++)
				result_weights[k].weight = bary_coord[k];
			return result_weights;
		}
	}

	//遍历所有边
	//思路：在所有重心为正的边中，找到距离最近的边, 并且采样点到该边的距离，一定小于到最近点的距离；
	T shortest_edge_dist = std::numeric_limits<T>::max();
	int closet_edge_index = -1;
	T edge_bray;
	for (int i = 0; i < _edges.size(); i++) {
		edge_bray = _edges[i].GetEdgeBaryCoord(sample_point);
		if (edge_bray > 0.0 && edge_bray < 1.0)
		{
			// 采样点 S ，线段两个端点 A B，垂点 P
			// P = (1-k) A + k B
			// SP * AB ==0
			// (Axaz-kAx+kBx -Sx) *(Bx - Ax )+(Ay-kAy+kBy -Sy) *(By - Ay ) = 0
			auto A = *_edges[i].v;
			auto B = *_edges[i].w;
			auto S = sample_point;
			auto AS = S - A;
			auto AB = B - A;
			auto AS2 = AS.norm2();
			auto ASAB = AS.dot(AB);

			T cos = ASAB / (AS.len() * AB.len());
			T sin = std::sqrt(1.0 - cos * cos);
			T point_edge_dist = AS.len() * sin;
			if (point_edge_dist < shortest_edge_dist)
			{
				shortest_edge_dist = point_edge_dist;
				closet_edge_index = i;
			}
		}
	}

	// 如果计算得到的线段的重心有负值，则查找离采样点最近的点，完全取该点权重为1
	double shortest_point_dist = std::numeric_limits<T>::max();
	dt::Vector2<T> closet_point;

	for (auto vertex : _vertices) {
		auto dist = vertex.dist(sample_point);
		if (dist < shortest_point_dist)
		{
			shortest_point_dist = dist;
			closet_point = vertex;
		}
	}

	// 最近点 OR 最近边
	if (shortest_edge_dist < shortest_point_dist)
	{
		std::cout << "Get Closet Edge:" << std::endl;
		edge_bray = _edges[closet_edge_index].GetEdgeBaryCoord(sample_point);

		for (int j = 0; j < _vertices.size(); j++) {
			if (almost_equal(*_edges[closet_edge_index].v, _vertices[j])) result_weights[0].vertex_index = j;
			if (almost_equal(*_edges[closet_edge_index].w, _vertices[j])) result_weights[1].vertex_index = j;
		}

		result_weights[0].weight = edge_bray;
		result_weights[1].weight = 1.0 - edge_bray;

		std::cout << "Closet Edge Index:" << closet_edge_index << std::endl;
		std::cout << "Closet Edge Bray Index:"<< closet_edge_index <<" Coord: (" << edge_bray << ", " << 1.0 - edge_bray << ")" << std::endl;
		std::cout << "Closet Edge Point V:(" << _edges[closet_edge_index].v->x << "," << _edges[closet_edge_index].v->y << " )" << std::endl;
		std::cout << "Closet Edge Point W:(" << _edges[closet_edge_index].w->x << "," << _edges[closet_edge_index].w->y << " )" << std::endl;
	}
	else
	{
		for (int j = 0; j < _vertices.size(); j++) {
			if (almost_equal(closet_point, _vertices[j]))
			{
				result_weights[0].vertex_index = j;
				result_weights[0].weight = 1.0;
			}
		}
	}

	return result_weights;
}

template<typename T>
const std::vector<typename Delaunay<T>::TriangleType>&
Delaunay<T>::getTriangles() const
{
	return _triangles;
}

template<typename T>
const std::vector<typename Delaunay<T>::EdgeType>&
Delaunay<T>::getEdges() const
{
	return _edges;
}

template<typename T>
const std::vector<typename Delaunay<T>::VertexType>&
Delaunay<T>::getVertices() const
{
	return _vertices;
}

template class Delaunay<float>;
template class Delaunay<double>;

} // namespace dt
