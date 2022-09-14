#include "edge.h"

namespace dt {

template<typename T>
Edge<T>::Edge(const VertexType &v1, const VertexType &v2) :
	v(&v1), w(&v2)
{}

template<typename T>
bool
Edge<T>::operator ==(const Edge<T> &e) const
{
	return (*(this->v) == *e.v && *(this->w) == *e.w) ||
			(*(this->v) == *e.w && *(this->w) == *e.v);
}

template<typename U>
std::ostream&
operator <<(std::ostream &str, const Edge<U> &e)
{
	return str << "Edge " << *e.v << ", " << *e.w;
}

template<typename T>
double Edge<T>::GetEdgeBaryCoord(const VertexType& sample_point) const
{
	auto A = *v;
	auto B = *w;
	auto S = sample_point;
	double k = (A.x * A.x + A.y * A.y + S.x * B.x - A.x * S.x + S.y * B.y - A.y * S.y - A.x * B.x - A.y * B.y) /
		(A.x * A.x + B.y * B.y + A.y * A.y + B.x * B.x - 2 * A.x * B.x - 2 * A.y * B.y);
	return 1.0-k;
};

template<typename T>
bool
Edge<T>::containsVertex(const VertexType& v_) const
{
	return almost_equal(*v, v_) || almost_equal(*w, v_);
}


template struct Edge<float>;
template struct Edge<double>;

} // namespace dt
