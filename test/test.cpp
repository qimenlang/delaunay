#define CATCH_CONFIG_MAIN
#include <random>
#include <catch2/catch_all.hpp>
#include "../delaunay/delaunay.h"

namespace dt {
	TEST_CASE("Test Edge get bary coord", "[Edge Test]") {

		std::vector<Vector2<double>> points;
		points.push_back(Vector2<double>{0.0, 0.0});
		points.push_back(Vector2<double>{1.0, 1.0});
		Vector2<double> Sample_point0 = { 0.0, 0.2 };
		Vector2<double> Sample_point1 = { 1.0, 0.0 };
		Vector2<double> Sample_point2 = { 0.0, 1.6 };

		Edge<double> edge = Edge<double>(points[0], points[1]);

		double Bary0 = edge.GetEdgeBaryCoord(Sample_point0);
		double Bary1 = edge.GetEdgeBaryCoord(Sample_point1);
		double Bary2 = edge.GetEdgeBaryCoord(Sample_point2);

		REQUIRE(almost_equal(Bary0, 0.9));
		REQUIRE(almost_equal(Bary1, 0.5));
		REQUIRE(almost_equal(Bary2, 0.2));
	}

	TEST_CASE("Test Vector2 minus", "[Vector2 Test]") {

		std::vector<Vector2<double>> points;
		points.push_back(Vector2<double>{0.0, 0.0});
		points.push_back(Vector2<double>{1.0, 0.0});
		points.push_back(Vector2<double>{0.0, 1.0});

		Vector2<double> point = {0.5,0.0};

		Triangle<double> tri = Triangle<double>(points[0],points[1],points[2]);
		auto BaryCentric = tri.GetBaryCentricCoord(point);
		for (auto coord : BaryCentric)
		{
			std::cout << coord << std::endl;
		}
		REQUIRE(BaryCentric[0] == 0.5);
		REQUIRE(BaryCentric[1] == 0.5);
		REQUIRE(BaryCentric[2] == 0.0);
	}

	TEST_CASE("Delaunay Get Bray Coord Weights", "[DelaunayTest]") {
		std::vector<Vector2<double>> points;
		points.push_back(Vector2<double>{0.0, 0.0});
		points.push_back(Vector2<double>{1.0, 0.0});
		points.push_back(Vector2<double>{0.0, 1.0});
		Delaunay<double> triangulation;
		const std::vector<Triangle<double>> triangles = triangulation.triangulate(points);

		dt::Vector2<double> sample_point = { 0.3,0.33 };
		dt::Vector2<double> lerp_point = { 0.0,0.0 };
		auto bray_coords = triangulation.getBrayCoord(sample_point);
		for (auto coord : bray_coords) {
			std::cout << coord.vertex_index << ":" << coord.weight << std::endl;
			lerp_point.x += points[coord.vertex_index].x * coord.weight;
			lerp_point.y += points[coord.vertex_index].y * coord.weight;
		}

		std::cout << "重心坐标插值点 :(" << lerp_point.x << "," << lerp_point.y << ")" << std::endl;
		REQUIRE(almost_equal(sample_point,lerp_point));
	}

	TEST_CASE("Delaunay triangulation should be able to triangulate 3 points as double", "[DelaunayTest]") {
		std::vector<Vector2<double>> points;
		points.push_back(Vector2<double>{0.0, 0.0});
		points.push_back(Vector2<double>{1.0, 0.0});
		points.push_back(Vector2<double>{0.0, 1.0});
		Delaunay<double> triangulation;
		const std::vector<Triangle<double>> triangles = triangulation.triangulate(points);
		REQUIRE(1 == triangles.size());
	}

	TEST_CASE("Delaunay triangulation should be able to handle duplicated 3 points as double", "[DelaunayTest]") {
		std::vector<Vector2<double>> points;
		points.push_back(Vector2<double>{0.0, 0.0});
		points.push_back(Vector2<double>{0.0, 0.0});
		points.push_back(Vector2<double>{1.0, 0.0});
		points.push_back(Vector2<double>{0.0, 1.0});
		Delaunay<double> triangulation;
		const std::vector<Triangle<double>> triangles = triangulation.triangulate(points);
		REQUIRE(1 == triangles.size());
		// const std::vector<Edge<double> > edges = triangulation.getEdges();
	}

	std::default_random_engine eng(std::random_device{}());

	TEST_CASE("Delaunay triangulation should be able to handle 10000 points as double", "[DelaunayTest]") {
		std::uniform_real_distribution<double> dist(0,
			std::numeric_limits<double>::max());
		constexpr size_t nb_pts = 1e4;
		std::vector<Vector2<double>> points(nb_pts);
		for (size_t i = 0; i < nb_pts; ++i)
		{
			points.at(i) = Vector2<double>{ dist(eng), dist(eng) };
		}
		REQUIRE(points.size() == nb_pts);
		Delaunay<double> triangulation;
		const std::vector<Triangle<double>> triangles = triangulation.triangulate(points);
	}

}
