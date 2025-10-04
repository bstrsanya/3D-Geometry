#ifndef TRIANGLES_HPP
#define TRIANGLES_HPP

#include <cmath>
#include <iostream>
#include <utility>

const double EPSILON = 1e-7;

// ------------------------------POINT_T---------------------------------------------

struct point_t
{
    double x_ = NAN;
    double y_ = NAN;
    double z_ = NAN;

    point_t () { };
    point_t (double x, double y, double z) : x_ { x }, y_ { y }, z_ { z } { };

    point_t operator+(const point_t& p) const { return { x_ + p.x_, y_ + p.y_, z_ + p.z_ }; }
    point_t operator-(const point_t& p) const { return { x_ - p.x_, y_ - p.y_, z_ - p.z_ }; }
    point_t operator*(double k) const { return { x_ * k, y_ * k, z_ * k }; }
    point_t operator/(double k) const { return { x_ / k, y_ / k, z_ / k };}
    void operator=(const point_t& p) { x_ = p.x_; y_ = p.y_; z_ = p.z_; }
    bool operator==(const point_t& p) const 
    { 
        return (std::fabs (x_ - p.x_) < EPSILON && 
                std::fabs (y_ - p.y_) < EPSILON &&
                std::fabs (z_ - p.z_) < EPSILON);
    }

    double get_x () const { return x_; }
    double get_y () const { return y_; }
    double get_z () const { return z_; }
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
    bool zero_vector () const 
    {
        return ((std::fabs (x_) < EPSILON) &&
                (std::fabs (y_) < EPSILON) && 
                (std::fabs (z_) < EPSILON));
    }
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
    point_t a_ {};
    point_t b_ {};
    point_t c_ {};

    vector_t N_; // the plane equation (N, X - a) = 0

    point_t p_min_ {};
    point_t p_max_ {};

public:
    triangle_t () { };
    triangle_t (const point_t& a, const point_t& b, const point_t& c);

    point_t get_a () const { return a_; }
    point_t get_b () const { return b_; }
    point_t get_c () const { return c_; }
    vector_t get_N () const { return N_; }

    double distance_point_plane_tr (const point_t& p) const;
    bool point_lie_in_plane_tr (const point_t& p) const;
    bool degenerate_tr () const;
    bool triangle_is_point () const { return ((a_ == b_) && (a_ == c_)); }
    bool triangle_is_line () const { return (degenerate_tr() && !triangle_is_point()); }
    bool triangle_lie_in_space (const point_t& p1, const point_t& p2) const;

    bool check_intersection (const triangle_t& other) const;
    bool check_same_sign_distance (const triangle_t& other) const;
    bool check_intersection_tr_of_line (const triangle_t& other) const ;
    std::pair<double, double> projection (char axis, const triangle_t& other) const;
    bool check_triangle_point (const point_t& p) const;
    bool check_triangle_line (const point_t& p1, const point_t& p2) const;
    static bool check_line_point (const point_t& line_p1, const point_t& line_p2, const point_t& p);
    static bool check_point_point (const point_t& p1, const point_t& p2);
    static bool check_line_line (const point_t& line1_p1, const point_t& line1_p2,
                      const point_t& line2_p1, const point_t& line2_p2);
    bool check_different_degeneracies (const triangle_t& other) const;
    static std::pair<point_t, point_t> select_ends_segment (const point_t& p1, 
                              const point_t& p2, const point_t& p3);
};

inline bool triangle_t::degenerate_tr () const
{
    return ((std::fabs(N_.get_x()) < EPSILON) &&
            (std::fabs(N_.get_y()) < EPSILON) && 
            (std::fabs(N_.get_z()) < EPSILON));
}

inline triangle_t::triangle_t (const point_t& a, const point_t& b, const point_t& c) : 
    a_(a), b_(b), c_(c), N_(vector_t ({ a, b }).cross_product ({ a, c }))
{
    p_min_.x_ = std::min(a_.x_, std::min (b_.x_, c_.x_));
    p_min_.y_ = std::min(a_.y_, std::min (b_.y_, c_.y_));
    p_min_.z_ = std::min(a_.z_, std::min (b_.z_, c_.z_));

    p_max_.x_ = std::max(a_.x_, std::max (b_.x_, c_.x_));
    p_max_.y_ = std::max(a_.y_, std::max (b_.y_, c_.y_));
    p_max_.z_ = std::max(a_.z_, std::max (b_.z_, c_.z_));
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
    return (std::fabs (distance_point_plane_tr (p)) < EPSILON);
}

inline bool triangle_t::triangle_lie_in_space(const point_t& p1, const point_t& p2) const
{
    double cube_min_x = std::min(p1.x_, p2.x_);
    double cube_max_x = std::max(p1.x_, p2.x_);

    double cube_min_y = std::min(p1.y_, p2.y_);
    double cube_max_y = std::max(p1.y_, p2.y_);

    double cube_min_z = std::min(p1.z_, p2.z_);
    double cube_max_z = std::max(p1.z_, p2.z_);

    bool overlap_x = !(p_max_.x_ < cube_min_x || p_min_.x_ > cube_max_x);
    bool overlap_y = !(p_max_.y_ < cube_min_y || p_min_.y_ > cube_max_y);
    bool overlap_z = !(p_max_.z_ < cube_min_z || p_min_.z_ > cube_max_z);

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
        std::pair<point_t, point_t> pair1 = select_ends_segment (a_, b_, c_);
        std::pair<point_t, point_t> pair2 = select_ends_segment (
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
        std::pair<point_t, point_t> pair = select_ends_segment (a_, b_, c_);
        return check_line_point (pair.first, pair.second, other.get_a ());
    }

    if (!degenerate_tr () && other.triangle_is_line ())
    {
        std::pair<point_t, point_t> pair = select_ends_segment (
                                        other.get_a (), other.get_b (), other.get_c ());
        return check_triangle_line (pair.first, pair.second);
    }

    if (!degenerate_tr () && other.triangle_is_point ())
    {
        return check_triangle_point (other.get_a ());
    }

    return false;
}

inline std::pair<point_t, point_t> triangle_t::select_ends_segment (const point_t& p1, 
                              const point_t& p2, const point_t& p3)
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

inline bool triangle_t::check_triangle_point (const point_t& p) const
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

inline bool triangle_t::check_triangle_line (const point_t& p1, const point_t& p2) const
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

    point_t p = p1 + ((p2 - p1) * t);

    return check_triangle_point (p);
}

inline bool triangle_t::check_line_line (const point_t& line1_p1, const point_t& line1_p2,
                      const point_t& line2_p1, const point_t& line2_p2)
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
  
    point_t line1_parameter1 = line1_p1 + (line1_p2 - line1_p1) * parameter_1;
    point_t line2_parameter2 = line2_p1 + (line2_p2 - line2_p1) * parameter_2;


    return ((vector_t{line1_parameter1, line2_parameter2}.zero_vector()) &&
            (parameter_1 >= -EPSILON) && (parameter_1 <= 1 + EPSILON)    &&
            (parameter_2 >= -EPSILON) && (parameter_2 <= 1 + EPSILON));
}

inline bool triangle_t::check_line_point (const point_t& line_p1, const point_t& line_p2, const point_t& p)
{
    vector_t vec1{ line_p1, line_p2 };
    vector_t vec2{ line_p1, p };

    if (!vec1.cross_product (vec2).zero_vector ())
    {
        return false;
    }

    return ((std::fmin (line_p1.get_x (), line_p2.get_x ()) - EPSILON <= p.get_x ()) &&
            (std::fmax (line_p1.get_x (), line_p2.get_x ()) + EPSILON >= p.get_x ()) &&

            (std::fmin (line_p1.get_y (), line_p2.get_y ()) - EPSILON <= p.get_y ()) &&
            (std::fmax (line_p1.get_y (), line_p2.get_y ()) + EPSILON >= p.get_y ()) &&

            (std::fmin (line_p1.get_z (), line_p2.get_z ()) - EPSILON <= p.get_z ()) &&
            (std::fmax (line_p1.get_z (), line_p2.get_z ()) + EPSILON >= p.get_z ()));
}

inline bool triangle_t::check_point_point (const point_t& p1, const point_t& p2)
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

                check_triangle_point (other.get_a ()) ||
                check_triangle_point (other.get_b ()) || 
                check_triangle_point (other.get_c ()) ||
                
                other.check_triangle_line (a_, b_) ||
                other.check_triangle_line (a_, c_) ||
                other.check_triangle_line (c_, b_) ||

                other.check_triangle_point (a_) ||
                other.check_triangle_point (b_) || 
                other.check_triangle_point (c_));
    }

    double D_x = std::fabs (D.get_x ());
    double D_y = std::fabs (D.get_y ());
    double D_z = std::fabs (D.get_z ());

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
    point_t point_1   = a_;
    point_t point_mid = b_;
    point_t point_2   = c_;

    // choose which axis to project on
    double project_point_1   = point_1.z_;
    double project_point_mid = point_mid.z_;
    double project_point_2   = point_2.z_;

    if (axis == 'x')
    {
        project_point_1   = point_1.x_;
        project_point_mid = point_mid.x_;
        project_point_2   = point_2.x_;
    }
    else if (axis == 'y')
    {
        project_point_1   = point_1.y_;
        project_point_mid = point_mid.y_;
        project_point_2   = point_2.y_;
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

std::istream& operator>> (std::istream& in, point_t& p);

std::istream& operator>> (std::istream& in, point_t& p)
{
    in >> p.x_ >> p.y_ >> p.z_;
    return in;
}

// ----------------------------------------------------------------------------------

#endif // TRIANGLES_HPP
