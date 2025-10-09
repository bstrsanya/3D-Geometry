#ifndef TRIANGLES_HPP
#define TRIANGLES_HPP

#include <cmath>
#include <iostream>
#include <utility>

const double EPSILON = 1e-7;

// ------------------------------VECTOR_T--------------------------------------------

class vector_t
{
public:
    double cor_x = NAN;
    double cor_y = NAN;
    double cor_z = NAN;

public:
    vector_t () = default;
    vector_t (double x, double y, double z) : cor_x { x }, cor_y { y }, cor_z { z } { };
    vector_t (const vector_t& a, const vector_t& b) : cor_x {b.cor_x - a.cor_x}, cor_y {b.cor_y - a.cor_y}, cor_z {b.cor_z - a.cor_z} {};

    vector_t operator+(const vector_t& p) const { return { cor_x + p.cor_x, cor_y + p.cor_y, cor_z + p.cor_z }; }
    vector_t operator-(const vector_t& p) const { return { cor_x - p.cor_x, cor_y - p.cor_y, cor_z - p.cor_z }; }
    vector_t operator*(double k)          const { return { cor_x * k, cor_y * k, cor_z * k }; }
    vector_t operator/(double k)          const { return { cor_x / k, cor_y / k, cor_z / k };}
    bool operator==(const vector_t& p) const;

    vector_t cross_product (const vector_t& b) const;
    double scalar_product (const vector_t& b) const;
    bool zero_vector () const;
};

inline vector_t vector_t::cross_product (const vector_t& b) const
{
    return vector_t { cor_y * b.cor_z - cor_z * b.cor_y,
                      cor_z * b.cor_x - cor_x * b.cor_z,
                      cor_x * b.cor_y - cor_y * b.cor_x};
}

inline double vector_t::scalar_product (const vector_t& b) const
{
    return cor_x * b.cor_x + cor_y * b.cor_y + cor_z * b.cor_z;
}

inline bool vector_t::operator==(const vector_t& p) const 
{
    return (std::fabs (cor_x - p.cor_x) < EPSILON && 
            std::fabs (cor_y - p.cor_y) < EPSILON &&
            std::fabs (cor_z - p.cor_z) < EPSILON);
}

inline bool vector_t::zero_vector () const 
{
    return ((std::fabs (cor_x) < EPSILON) &&
            (std::fabs (cor_y) < EPSILON) && 
            (std::fabs (cor_z) < EPSILON));
}

// ----------------------------------------------------------------------------------

// ------------------------------TRIANGLE_T------------------------------------------

class triangle_t
{
    vector_t a_ {};
    vector_t b_ {};
    vector_t c_ {};

    vector_t N_; // the plane equation (N, X - a) = 0

    vector_t p_min_ {};
    vector_t p_max_ {};

public:
    triangle_t () { };
    triangle_t (const vector_t& a, const vector_t& b, const vector_t& c);

    vector_t get_a () const { return a_; }
    vector_t get_b () const { return b_; }
    vector_t get_c () const { return c_; }
    vector_t get_N () const { return N_; }

    double distance_point_plane_tr (const vector_t& p) const;
    bool point_lie_in_plane_tr (const vector_t& p) const;
    bool degenerate_tr () const;
    bool triangle_is_point () const { return ((a_ == b_) && (a_ == c_)); }
    bool triangle_is_line () const { return (degenerate_tr() && !triangle_is_point()); }
    bool triangle_lie_in_space (const vector_t& p1, const vector_t& p2) const;

    bool check_intersection (const triangle_t& other) const;
    bool check_same_sign_distance (const triangle_t& other) const;
    bool check_intersection_tr_of_line (const triangle_t& other) const ;
    std::pair<double, double> projection (char axis, const triangle_t& other) const;
    bool check_triangle_point (const vector_t& p) const;
    bool check_triangle_line (const vector_t& p1, const vector_t& p2) const;
    static bool check_line_point (const vector_t& line_p1, const vector_t& line_p2, const vector_t& p);
    static bool check_point_point (const vector_t& p1, const vector_t& p2);
    static bool check_line_line (const vector_t& line1_p1, const vector_t& line1_p2,
                      const vector_t& line2_p1, const vector_t& line2_p2);
    bool check_different_degeneracies (const triangle_t& other) const;
    static std::pair<vector_t, vector_t> select_ends_segment (const vector_t& p1, 
                              const vector_t& p2, const vector_t& p3);
};

inline bool triangle_t::degenerate_tr () const
{
    return ((std::fabs(N_.cor_x) < EPSILON) &&
            (std::fabs(N_.cor_y) < EPSILON) && 
            (std::fabs(N_.cor_z) < EPSILON));
}

inline triangle_t::triangle_t (const vector_t& a, const vector_t& b, const vector_t& c) : 
    a_(a), b_(b), c_(c), N_(vector_t ({ a, b }).cross_product ({ a, c }))
{
    p_min_.cor_x = std::min(a_.cor_x, std::min (b_.cor_x, c_.cor_x));
    p_min_.cor_y = std::min(a_.cor_y, std::min (b_.cor_y, c_.cor_y));
    p_min_.cor_z = std::min(a_.cor_z, std::min (b_.cor_z, c_.cor_z));

    p_max_.cor_x = std::max(a_.cor_x, std::max (b_.cor_x, c_.cor_x));
    p_max_.cor_y = std::max(a_.cor_y, std::max (b_.cor_y, c_.cor_y));
    p_max_.cor_z = std::max(a_.cor_z, std::max (b_.cor_z, c_.cor_z));
}

inline double triangle_t::distance_point_plane_tr (const vector_t& p) const
{
    vector_t vec { a_, p };
    double res  = N_.scalar_product (vec);
    double norm = std::sqrt (N_.scalar_product (N_));
    return res / norm;
}

inline bool triangle_t::point_lie_in_plane_tr (const vector_t& p) const
{
    return (std::fabs (distance_point_plane_tr (p)) < EPSILON);
}

inline bool triangle_t::triangle_lie_in_space(const vector_t& p1, const vector_t& p2) const
{
    double cube_min_x = std::min(p1.cor_x, p2.cor_x);
    double cube_max_x = std::max(p1.cor_x, p2.cor_x);

    double cube_min_y = std::min(p1.cor_y, p2.cor_y);
    double cube_max_y = std::max(p1.cor_y, p2.cor_y);

    double cube_min_z = std::min(p1.cor_z, p2.cor_z);
    double cube_max_z = std::max(p1.cor_z, p2.cor_z);

    bool overlap_x = !(p_max_.cor_x < cube_min_x || p_min_.cor_x > cube_max_x);
    bool overlap_y = !(p_max_.cor_y < cube_min_y || p_min_.cor_y > cube_max_y);
    bool overlap_z = !(p_max_.cor_z < cube_min_z || p_min_.cor_z > cube_max_z);

    return overlap_x && overlap_y && overlap_z;
}

inline bool triangle_t::check_intersection (const triangle_t& other) const
{
    if (other.check_same_sign_distance (*this) || 
        check_same_sign_distance (other))
    {
        return false; // one of the triangles lies completely in the half-plane of
                      // the other
    }

    // non-degenerate triangles
    if (!degenerate_tr () && !other.degenerate_tr ())
        return check_intersection_tr_of_line (other);

    // the same kind of degeneracy
    if (triangle_is_line () && other.triangle_is_line ())
    {
        std::pair<vector_t, vector_t> pair1 = select_ends_segment (a_, b_, c_);
        std::pair<vector_t, vector_t> pair2 = select_ends_segment (
                                        other.get_a (), other.get_b (), other.get_c ());
        return check_line_line (pair1.first, pair1.second, pair2.first, pair2.second);
    }

    if (triangle_is_point () && other.triangle_is_point ())
        return check_point_point (a_, other.get_a ());

    // different kinds of degeneracy
    return (check_different_degeneracies (other) || 
            other.check_different_degeneracies (*this));
}

inline bool triangle_t::check_different_degeneracies (const triangle_t& other) const 
{
    if (triangle_is_line () && other.triangle_is_point ())
    {
        std::pair<vector_t, vector_t> pair = select_ends_segment (a_, b_, c_);
        return check_line_point (pair.first, pair.second, other.get_a ());
    }

    if (!degenerate_tr () && other.triangle_is_line ())
    {
        std::pair<vector_t, vector_t> pair = select_ends_segment (
                                        other.get_a (), other.get_b (), other.get_c ());
        return check_triangle_line (pair.first, pair.second);
    }

    if (!degenerate_tr () && other.triangle_is_point ())
    {
        return check_triangle_point (other.get_a ());
    }

    return false;
}

inline std::pair<vector_t, vector_t> triangle_t::select_ends_segment (const vector_t& p1, 
                              const vector_t& p2, const vector_t& p3)
{
    // c_ lies between a_ & b_
    if (check_line_point (p1, p2, p3)) 
        return {p1, p2};

    // b_ lies between a_ & c_
    else if (check_line_point (p1, p3, p2))
        return {p1, p3};

    // a_ lies between b_ & c_
    return {p2, p3};
};

inline bool triangle_t::check_triangle_point (const vector_t& p) const
{
    // guaranteed that point already lies in plane of triangle
    double res_1 = (vector_t{a_, b_}.cross_product(
                    vector_t{a_, p})).scalar_product(N_); 

    double res_2 = (vector_t{b_, c_}.cross_product(
                    vector_t{b_, p})).scalar_product(N_); 
    
    double res_3 = (vector_t{c_, a_}.cross_product(
                    vector_t{c_, p})).scalar_product(N_); 

    return ((res_1 >= -EPSILON && res_2 >= -EPSILON && res_3 >= -EPSILON) || 
            (res_1 <= EPSILON && res_2 <= EPSILON && res_3 <= EPSILON));
}

inline bool triangle_t::check_triangle_line (const vector_t& p1, const vector_t& p2) const
{
    // guaranteed that line intersects plane of triangle (or lies)

    // lies
    if (N_.scalar_product(vector_t{p1, p2}) == 0)
    {
         return (check_line_line (p1, p2, a_, b_) ||
                 check_line_line (p1, p2, b_, c_) ||
                 check_line_line (p1, p2, c_, a_) ||
                 check_triangle_point (p1) || check_triangle_point (p2));
    }
    
    // intersects
    double t = - N_.scalar_product(vector_t{a_, p1}) /
                 N_.scalar_product(vector_t{p1, p2});

    vector_t p = p1 + ((p2 - p1) * t);

    return check_triangle_point (p);
}

inline bool triangle_t::check_line_line (const vector_t& line1_p1, const vector_t& line1_p2,
                      const vector_t& line2_p1, const vector_t& line2_p2)
{
    vector_t u{ line1_p1, line1_p2 }; // guiding vector line 1
    vector_t v{ line2_p1, line2_p2 }; // guiding vector line 2
    vector_t w{ line1_p1, line2_p1 }; // vector between line 1 & line 2

    // parallel lines
    if ((u.cross_product(v)).zero_vector()) 
    {
        return (check_line_point(line1_p1, line1_p2, line2_p1) || 
                check_line_point(line1_p1, line1_p2, line2_p2) ||
                check_line_point(line2_p1, line2_p2, line1_p1) || 
                check_line_point(line2_p1, line2_p2, line1_p2));
    }

    double u_u = u.scalar_product(u);
    double u_v = u.scalar_product(v);
    double v_v = v.scalar_product(v);
    double u_w = u.scalar_product(w);
    double v_w = v.scalar_product(w);

    double denominator = u_u * v_v - u_v * u_v;
    double parameter_1 = (v_v * u_w - u_v * v_w) / denominator;
    double parameter_2 = (u_v * u_w - u_u * v_w) / denominator;
  
    vector_t line1_parameter1 = line1_p1 + (line1_p2 - line1_p1) * parameter_1;
    vector_t line2_parameter2 = line2_p1 + (line2_p2 - line2_p1) * parameter_2;


    return ((vector_t{line1_parameter1, line2_parameter2}.zero_vector()) &&
            (parameter_1 >= -EPSILON) && (parameter_1 <= 1 + EPSILON)    &&
            (parameter_2 >= -EPSILON) && (parameter_2 <= 1 + EPSILON));
}

inline bool triangle_t::check_line_point (const vector_t& line_p1, const vector_t& line_p2, const vector_t& p)
{
    vector_t vec1{ line_p1, line_p2 };
    vector_t vec2{ line_p1, p };

    if (!vec1.cross_product (vec2).zero_vector ())
    {
        return false;
    }

    return ((std::fmin (line_p1.cor_x, line_p2.cor_x) - EPSILON <= p.cor_x) &&
            (std::fmax (line_p1.cor_x, line_p2.cor_x) + EPSILON >= p.cor_x) &&

            (std::fmin (line_p1.cor_y, line_p2.cor_y) - EPSILON <= p.cor_y) &&
            (std::fmax (line_p1.cor_y, line_p2.cor_y) + EPSILON >= p.cor_y) &&

            (std::fmin (line_p1.cor_z, line_p2.cor_z) - EPSILON <= p.cor_z) &&
            (std::fmax (line_p1.cor_z, line_p2.cor_z) + EPSILON >= p.cor_z));
}

inline bool triangle_t::check_point_point (const vector_t& p1, const vector_t& p2)
{
    return (p1 == p2);
}

inline bool triangle_t::check_same_sign_distance (const triangle_t& other) const
{
    if (other.degenerate_tr())
    {
        return false;
    }
    double distance_1 = other.distance_point_plane_tr (a_);
    double distance_2 = other.distance_point_plane_tr (b_);
    double distance_3 = other.distance_point_plane_tr (c_);

    return ((distance_1 > EPSILON && distance_2 > EPSILON && distance_3 > EPSILON) || 
            (distance_1 < -EPSILON && distance_2 < -EPSILON && distance_3 < -EPSILON));
}

inline bool triangle_t::check_intersection_tr_of_line (const triangle_t& other) const
{
    // D = direction of the common line
    vector_t D = N_.cross_product (other.get_N ());

    // triangles lie in same plane
    if (D.zero_vector ())
    {
        return (check_triangle_line (other.get_a (), other.get_b ()) ||
                check_triangle_line (other.get_a (), other.get_c ()) ||
                check_triangle_line (other.get_c (), other.get_b ()) ||

                other.check_triangle_line (a_, b_) ||
                other.check_triangle_line (a_, c_) ||
                other.check_triangle_line (c_, b_));
    }

    double D_x = std::fabs (D.cor_x);
    double D_y = std::fabs (D.cor_y);
    double D_z = std::fabs (D.cor_z);

    char axis = 'z';

    if ((D_y - D_x) < EPSILON && (D_z - D_x) < EPSILON)
        axis = 'x';
    else if ((D_x - D_y) < EPSILON && (D_z - D_y) < EPSILON)
        axis = 'y';

    std::pair<double, double> pair_1 = projection (axis, other);
    double t1 = pair_1.first;
    double t2 = pair_1.second;

    std::pair<double, double> pair_2 = other.projection (axis, *this);
    double t3 = pair_2.first;
    double t4 = pair_2.second;

    return ((std::min (t1, t2) <= std::max (t3, t4)) && 
            (std::min (t3, t4) <= std::max (t1, t2)));
}

inline std::pair<double, double> triangle_t::projection (char axis, const triangle_t& other) const
{
    // point_1 and point_2 lie on the same side of tr_2 and that point_mid lies on
    // the other side
    vector_t point_1   = a_;
    vector_t point_mid = b_;
    vector_t point_2   = c_;

    // choose which axis to project on
    double project_point_1   = point_1.cor_z;
    double project_point_mid = point_mid.cor_z;
    double project_point_2   = point_2.cor_z;

    if (axis == 'x')
    {
        project_point_1   = point_1.cor_x;
        project_point_mid = point_mid.cor_x;
        project_point_2   = point_2.cor_x;
    }
    else if (axis == 'y')
    {
        project_point_1   = point_1.cor_y;
        project_point_mid = point_mid.cor_y;
        project_point_2   = point_2.cor_y;
    }

    // if two points lie on a plane tr_2
    if (other.point_lie_in_plane_tr (a_) && 
        other.point_lie_in_plane_tr (b_))
    {
        return { project_point_1, project_point_mid };
    }

    if (other.point_lie_in_plane_tr (b_) && 
        other.point_lie_in_plane_tr (c_))
    {
        return { project_point_2, project_point_mid };
    }
    if (other.point_lie_in_plane_tr (a_) && 
        other.point_lie_in_plane_tr (c_))
    {
        return { project_point_1, project_point_2 };
    }

    // By default, the middle point is point B, below are 3 cases where point B cannot be middle

    if (other.distance_point_plane_tr (b_) * 
        other.distance_point_plane_tr (c_) > 0)
    {
        std::swap (project_point_1, project_point_mid);
        std::swap (point_1, point_mid);
    }
    else if (other.distance_point_plane_tr (a_) * 
             other.distance_point_plane_tr (b_) > 0)
    {
        std::swap (project_point_2, project_point_mid);
        std::swap (point_2, point_mid);
    }

    if (other.point_lie_in_plane_tr (b_) && 
        other.distance_point_plane_tr (a_) * 
        other.distance_point_plane_tr (c_) < 0)
    {
        std::swap (project_point_1, project_point_mid);
        std::swap (point_1, point_mid);
    }

    // counting projections on the selected axis
    double t1 = project_point_1 + 
                (project_point_mid - project_point_1) * (other.distance_point_plane_tr (point_1) / 
                (other.distance_point_plane_tr (point_1) - other.distance_point_plane_tr (point_mid)));

    double t2 = project_point_2 + 
                (project_point_mid - project_point_2) * (other.distance_point_plane_tr (point_2) / 
                (other.distance_point_plane_tr (point_2) - other.distance_point_plane_tr (point_mid)));

    return { t1, t2 };
}

// ----------------------------------------------------------------------------------

// ------------------------------OTHER_FUNC------------------------------------------

std::istream& operator>> (std::istream& in, vector_t& p);

std::istream& operator>> (std::istream& in, vector_t& p)
{
    in >> p.cor_x >> p.cor_y >> p.cor_z;
    return in;
}

// ----------------------------------------------------------------------------------

#endif // TRIANGLES_HPP
