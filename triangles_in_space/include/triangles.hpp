#ifndef TRIANGLES_HPP
#define TRIANGLES_HPP

#include <cmath>
#include <iostream>
#include <utility>

const double epsilon                = 1e-7;
const double epsilon_for_degenerate = 1e-32;

// ------------------------------POINT_T---------------------------------------------

struct point_t
{
    double x_ = NAN;
    double y_ = NAN;
    double z_ = NAN;

    point_t () { };
    point_t (double x, double y, double z) : x_ { x }, y_ { y }, z_ { z } { };
};

// ----------------------------------------------------------------------------------

// ------------------------------VECTOR_T--------------------------------------------

class vector_t
{
    double x_ = NAN;
    double y_ = NAN;
    double z_ = NAN;

public:
    vector_t () { };
    vector_t (double x, double y, double z) : x_ { x }, y_ { y }, z_ { z } { };
    vector_t (const point_t& a, const point_t& b);

    double get_x () const { return x_; }
    double get_y () const { return y_; }
    double get_z () const { return z_; }

    vector_t cross_product (const vector_t& b) const;
    double scalar_product (const vector_t& b) const;
};

inline vector_t::vector_t (const point_t& a, const point_t& b)
{
    x_ = b.x_ - a.x_;
    y_ = b.y_ - a.y_;
    z_ = b.z_ - a.z_;
}

inline vector_t vector_t::cross_product (const vector_t& b) const
{
    vector_t c;
    c.x_ = y_ * b.z_ - z_ * b.y_;
    c.y_ = z_ * b.x_ - x_ * b.z_;
    c.z_ = x_ * b.y_ - y_ * b.x_;
    return c;
}

inline double vector_t::scalar_product (const vector_t& b) const
{
    return x_ * b.x_ + y_ * b.y_ + z_ * b.z_;
}

// ----------------------------------------------------------------------------------

// ------------------------------TRIANGLE_T------------------------------------------

class triangle_t
{
    point_t a_;
    point_t b_;
    point_t c_;

    vector_t N_; // the plane equation (N, X - a) = 0

public:
    triangle_t () { };
    triangle_t (const point_t& a, const point_t& b, const point_t& c);

    point_t get_a () const { return a_; }
    point_t get_b () const { return b_; }
    point_t get_c () const { return c_; }
    vector_t get_N () const { return N_; }

    double distance_point_plane_tr (const point_t& p) const;
    bool point_lie_in_plane_tr (const point_t& p) const;
};

inline triangle_t::triangle_t (const point_t& a, const point_t& b, const point_t& c)
{
    vector_t vec_1 { a, b };
    vector_t vec_2 { a, c };
    vector_t cross = vec_1.cross_product (vec_2);
    double res     = cross.scalar_product (cross);
    if (res < epsilon_for_degenerate)
    {
        std::cerr << "Degenerate triangles are not accepted!\n";
        exit (1);
    }
    a_ = a;
    b_ = b;
    c_ = c;
    N_ = vector_t ({ a, b }).cross_product ({ a, c });
}

inline double triangle_t::distance_point_plane_tr (const point_t& p) const
{
    vector_t vec { a_, p };
    double res  = N_.scalar_product (vec);
    double norm = std::sqrt (N_.scalar_product (N_));
    return res / norm;
}

inline bool triangle_t::point_lie_in_plane_tr (const point_t& p) const
{
    return (std::fabs (distance_point_plane_tr (p)) < epsilon);
}

// ----------------------------------------------------------------------------------

// ------------------------------DECLARATIONS_FUNC-----------------------------------

std::istream& operator>> (std::istream& in, point_t& p);

bool check_intersection (const triangle_t& tr_1, const triangle_t& tr_2);
bool check_same_sign_distance (const triangle_t& tr_1, const triangle_t& tr_2);
bool check_intersection_tr_of_line (const triangle_t& tr_1, const triangle_t& tr_2);
std::pair<double, double> projection (char axis, const triangle_t& tr_1, const triangle_t& tr_2);

// ----------------------------------------------------------------------------------

#endif // TRIANGLES_HPP
