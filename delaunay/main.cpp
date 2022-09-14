#include <iostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <array>
#include <random>
#include <chrono>
#include <limits>

#include <SFML/Graphics.hpp>

#include "vector2.h"
#include "triangle.h"
#include "delaunay.h"
int main(int argc, char * argv[])
{
	int numberPoints = 40;
	if (argc>1)
	{
		numberPoints = atoi(argv[1]);
	}

	std::default_random_engine eng(std::random_device{}());
	std::uniform_real_distribution<double> dist_w(0, 800);
	std::uniform_real_distribution<double> dist_h(0, 600);

	std::cout << "Generating " << numberPoints << " random points" << std::endl;

	std::vector<dt::Vector2<double>> points;
	for(int i = 0; i < numberPoints; ++i) {
		points.push_back(dt::Vector2<double>{dist_w(eng), dist_h(eng)});
	}
	//points.push_back(dt::Vector2<double>{40.0,70.0});
	//points.push_back(dt::Vector2<double>{70.0, 40.0});
	//points.push_back(dt::Vector2<double>{57.0, 57.0});

	dt::Delaunay<double> triangulation;
	const auto delaunay_start = std::chrono::high_resolution_clock::now();
	const std::vector<dt::Triangle<double>> triangles = triangulation.triangulate(points);
	const auto delaunay_end = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double> delaunay_cost = delaunay_end - delaunay_start;

	std::cout << triangles.size() << " triangles generated in " << delaunay_cost.count()
			<< "s\n";

	dt::Vector2<double> SamplePoint{ 50.0,50.0};
	//dt::Vector2<double> SamplePoint{ dist_w(eng), dist_h(eng) };

	std::cout << "Sample Point: "<<SamplePoint.x<<","<<SamplePoint.y << std::endl;

	std::cout << " Located In trangle blew :" << std::endl;

	//判断采样点位于哪个三角形内部
	std::array<double, 3> bary_coord;

	for (int i = 0;i<triangles.size();i++)
	{
		bary_coord = triangles[i].GetBaryCentricCoord(SamplePoint);
		if(bary_coord[0]>0&& bary_coord[1] > 0&& bary_coord[2] > 0)
			std::cout <<i<<":(" <<bary_coord[0] << "," << bary_coord[1] << "," << bary_coord[2]<<")" << std::endl;
	}

	auto bray_coords =  triangulation.getBrayCoord(SamplePoint);

	//使用重心坐标插值得到新的点，与采样点相同则正确
	dt::Vector2<double> lerp_point = {0.0,0.0};
	for (auto coord : bray_coords) {
		if (coord.weight <= 0) break;
		std::cout << coord.vertex_index << ":" << coord.weight << ":("<< points[coord.vertex_index].x<<","<< points[coord.vertex_index].y<<")"<<std::endl;
		lerp_point.x += points[coord.vertex_index].x * coord.weight;
		lerp_point.y += points[coord.vertex_index].y * coord.weight;
	}
	std::cout << "重心坐标插值点 :(" << lerp_point.x << "," << lerp_point.y << ")" << std::endl;

	//如果采样点不在任何三角形内，则判断距离哪条边最近，并计算采样点在该边上的重心
	
	const auto bray_end = std::chrono::high_resolution_clock::now();
	const std::chrono::duration<double> bray_cost = bray_end - delaunay_end;
	std::cout << triangles.size() << " compute Bray Coord Cost: " << bray_cost.count()
		<< "s\n";

	const std::vector<dt::Edge<double>> edges = triangulation.getEdges();

	// SFML window
	sf::RenderWindow window(sf::VideoMode(800, 600), "Delaunay triangulation");
	window.setFramerateLimit(1);

	// Draw The Sample Point
	sf::RectangleShape s{ sf::Vector2f(4, 4) };
	s.setPosition(static_cast<float>(SamplePoint.x),static_cast<float>(SamplePoint.y));
	window.draw(s);

	// Transform each points of each vector as a rectangle
	for(const auto p : points) {
		sf::RectangleShape s{sf::Vector2f(4, 4)};
		s.setPosition(static_cast<float>(p.x), static_cast<float>(p.y));
		window.draw(s);
	}
	

	std::vector<std::array<sf::Vertex, 2> > lines;
	for(const auto &e : edges) {
		const std::array<sf::Vertex, 2> line{{
			sf::Vertex(sf::Vector2f(
					static_cast<float>(e.v->x + 2.),
					static_cast<float>(e.v->y + 2.))),
			sf::Vertex(sf::Vector2f(
					static_cast<float>(e.w->x + 2.),
					static_cast<float>(e.w->y + 2.))),
		}};
		window.draw(line.data(), 2, sf::Lines);
	}

	window.display();

	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
		}
	}

	return 0;
}
