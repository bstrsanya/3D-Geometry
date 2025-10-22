#ifndef TRIANGLES_HPP
#define TRIANGLES_HPP

#include <cmath>
#include <iostream>
#include <utility>

// ------------------------------VECTOR_T--------------------------------------------

namespace Triangle {

template <typename TypeNum>
const TypeNum EPSILON = 1e-7;

template <typename TypeNum>
class vector_t
{
public:
    TypeNum cor_x = NAN;
    TypeNum cor_y = NAN;
    TypeNum cor_z = NAN;

public:
    vector_t () = default;
    vector_t (TypeNum x, TypeNum y, TypeNum z) : cor_x { x }, cor_y { y }, cor_z { z } { };
    vector_t (const vector_t<TypeNum>& a, const vector_t<TypeNum>& b) 
            : cor_x {b.cor_x - a.cor_x}, cor_y {b.cor_y - a.cor_y}, cor_z {b.cor_z - a.cor_z} {};

    vector_t operator+(const vector_t<TypeNum>& p) const 
                                { return { cor_x + p.cor_x, cor_y + p.cor_y, cor_z + p.cor_z }; }
    vector_t operator-(const vector_t<TypeNum>& p) const 
                                { return { cor_x - p.cor_x, cor_y - p.cor_y, cor_z - p.cor_z }; }
    vector_t operator*(TypeNum k) const { return { cor_x * k, cor_y * k, cor_z * k }; }
    vector_t operator/(TypeNum k) const { return { cor_x / k, cor_y / k, cor_z / k };}
    bool operator==(const vector_t& p) const;

    vector_t<TypeNum> cross_product (const vector_t<TypeNum>& b) const;
    TypeNum scalar_product (const vector_t<TypeNum>& b) const;
    bool zero_vector () const;
};

template <typename TypeNum>
inline vector_t<TypeNum> vector_t<TypeNum>::cross_product (const vector_t<TypeNum>& b) const
{
    return vector_t { cor_y * b.cor_z - cor_z * b.cor_y,
                      cor_z * b.cor_x - cor_x * b.cor_z,
                      cor_x * b.cor_y - cor_y * b.cor_x };
}

template <typename TypeNum>
inline TypeNum vector_t<TypeNum>::scalar_product (const vector_t<TypeNum>& b) const
{
    return cor_x * b.cor_x + cor_y * b.cor_y + cor_z * b.cor_z;
}

template <typename TypeNum>
inline bool vector_t<TypeNum>::operator==(const vector_t<TypeNum>& p) const 
{
    return (std::fabs (cor_x - p.cor_x) < EPSILON<TypeNum> && 
            std::fabs (cor_y - p.cor_y) < EPSILON<TypeNum> &&
            std::fabs (cor_z - p.cor_z) < EPSILON<TypeNum>);
}

template <typename TypeNum>
inline bool vector_t<TypeNum>::zero_vector () const 
{
    return ((std::fabs (cor_x) < EPSILON<TypeNum>) &&
            (std::fabs (cor_y) < EPSILON<TypeNum>) && 
            (std::fabs (cor_z) < EPSILON<TypeNum>));
}
} // ClassVector

// ----------------------------------------------------------------------------------

// ------------------------------TRIANGLE_T------------------------------------------

namespace Triangle {

template <typename TypeNum>
class triangle_t
{
    vector_t<TypeNum> a_ {};
    vector_t<TypeNum> b_ {};
    vector_t<TypeNum> c_ {};

    vector_t<TypeNum> N_; // the plane equation (N, X - a) = 0

    vector_t<TypeNum> p_min_ {};
    vector_t<TypeNum> p_max_ {};

    enum class axis_t
    {
        X_axis,
        Y_axis,
        Z_axis
    };

public:
    triangle_t () = default;
    triangle_t (const vector_t<TypeNum>& a, const vector_t<TypeNum>& b, const vector_t<TypeNum>& c);

    const vector_t<TypeNum>& get_a () const { return a_; }
    const vector_t<TypeNum>& get_b () const { return b_; }
    const vector_t<TypeNum>& get_c () const { return c_; }
    const vector_t<TypeNum>& get_N () const { return N_; }

    TypeNum distance_point_plane_tr (const vector_t<TypeNum>& p) const;
    bool point_lie_in_plane_tr (const vector_t<TypeNum>& p) const;
    bool degenerate_tr () const;
    bool triangle_lie_in_space (const vector_t<TypeNum>& p1, const vector_t<TypeNum>& p2) const;

    bool check_intersection_tr_tr (const triangle_t& other) const;
    bool check_intersection_tr_point (const vector_t<TypeNum>& p) const;
    bool check_intersection_tr_line (const vector_t<TypeNum>& p1, 
                                     const vector_t<TypeNum>& p2) const;
    static bool check_intersection_line_point (const vector_t<TypeNum>& line_p1, 
                                               const vector_t<TypeNum>& line_p2, 
                                               const vector_t<TypeNum>& p);
    static bool check_intersection_point_point (const vector_t<TypeNum>& p1, 
                                                const vector_t<TypeNum>& p2);
    static bool check_intersection_line_line (const vector_t<TypeNum>& line1_p1, 
                                              const vector_t<TypeNum>& line1_p2,
                                              const vector_t<TypeNum>& line2_p1, 
                                              const vector_t<TypeNum>& line2_p2);
    
private:
    static std::pair<vector_t<TypeNum>, vector_t<TypeNum>> select_ends_segment (
                                                              const vector_t<TypeNum>& p1, 
                                                              const vector_t<TypeNum>& p2, 
                                                              const vector_t<TypeNum>& p3);
    std::pair<TypeNum, TypeNum> projection (axis_t axis, const triangle_t& other) const;
    bool check_different_degeneracies (const triangle_t& other) const;
    bool triangle_is_point () const { return ((a_ == b_) && (a_ == c_)); }
    bool triangle_is_line () const { return (degenerate_tr() && !triangle_is_point()); }
    bool check_same_sign_distance (const triangle_t& other) const;
    bool check_intersection_tr_of_line (const triangle_t& other) const;
};

template <typename TypeNum>
inline bool triangle_t<TypeNum>::degenerate_tr () const
{
    return ((std::fabs(N_.cor_x) < EPSILON<TypeNum>) &&
            (std::fabs(N_.cor_y) < EPSILON<TypeNum>) && 
            (std::fabs(N_.cor_z) < EPSILON<TypeNum>));
}

template <typename TypeNum>
inline triangle_t<TypeNum>::triangle_t (const vector_t<TypeNum>& a, 
            const vector_t<TypeNum>& b, const vector_t<TypeNum>& c) : 
    a_(a), b_(b), c_(c), N_(vector_t<TypeNum> ({ a, b }).cross_product ({ a, c }))
{
    p_min_.cor_x = std::min(a_.cor_x, std::min (b_.cor_x, c_.cor_x));
    p_min_.cor_y = std::min(a_.cor_y, std::min (b_.cor_y, c_.cor_y));
    p_min_.cor_z = std::min(a_.cor_z, std::min (b_.cor_z, c_.cor_z));

    p_max_.cor_x = std::max(a_.cor_x, std::max (b_.cor_x, c_.cor_x));
    p_max_.cor_y = std::max(a_.cor_y, std::max (b_.cor_y, c_.cor_y));
    p_max_.cor_z = std::max(a_.cor_z, std::max (b_.cor_z, c_.cor_z));
}

template <typename TypeNum>
inline TypeNum triangle_t<TypeNum>::distance_point_plane_tr (const vector_t<TypeNum>& p) const
{
    vector_t vec { a_, p };
    TypeNum res  = N_.scalar_product (vec);
    TypeNum norm = std::sqrt (N_.scalar_product (N_));
    return res / norm;
}

template <typename TypeNum>
inline bool triangle_t<TypeNum>::point_lie_in_plane_tr (const vector_t<TypeNum>& p) const
{
    return (std::fabs (distance_point_plane_tr (p)) < EPSILON<TypeNum>);
}

template <typename TypeNum>
inline bool triangle_t<TypeNum>::triangle_lie_in_space(const vector_t<TypeNum>& p1, 
                                                       const vector_t<TypeNum>& p2) const
{
    TypeNum cube_min_x = std::min(p1.cor_x, p2.cor_x);
    TypeNum cube_max_x = std::max(p1.cor_x, p2.cor_x);

    TypeNum cube_min_y = std::min(p1.cor_y, p2.cor_y);
    TypeNum cube_max_y = std::max(p1.cor_y, p2.cor_y);

    TypeNum cube_min_z = std::min(p1.cor_z, p2.cor_z);
    TypeNum cube_max_z = std::max(p1.cor_z, p2.cor_z);

    bool overlap_x = !(p_max_.cor_x < cube_min_x || p_min_.cor_x > cube_max_x);
    bool overlap_y = !(p_max_.cor_y < cube_min_y || p_min_.cor_y > cube_max_y);
    bool overlap_z = !(p_max_.cor_z < cube_min_z || p_min_.cor_z > cube_max_z);

    return overlap_x && overlap_y && overlap_z;
}

template <typename TypeNum>
inline bool triangle_t<TypeNum>::check_intersection_tr_tr (const triangle_t& other) const
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
        std::pair<vector_t<TypeNum>, vector_t<TypeNum>> pair1 = select_ends_segment (a_, b_, c_);
        std::pair<vector_t<TypeNum>, vector_t<TypeNum>> pair2 = select_ends_segment (
                                        other.get_a (), other.get_b (), other.get_c ());
        return check_intersection_line_line (pair1.first, pair1.second, pair2.first, pair2.second);
    }

    if (triangle_is_point () && other.triangle_is_point ())
        return check_intersection_point_point (a_, other.get_a ());

    // different kinds of degeneracy
    return (check_different_degeneracies (other) || 
            other.check_different_degeneracies (*this));
}

template <typename TypeNum>
inline bool triangle_t<TypeNum>::check_different_degeneracies (const triangle_t& other) const 
{
    if (triangle_is_line () && other.triangle_is_point ())
    {
        std::pair<vector_t<TypeNum>, vector_t<TypeNum>> pair = select_ends_segment (a_, b_, c_);
        return check_intersection_line_point (pair.first, pair.second, other.get_a ());
    }

    if (!degenerate_tr () && other.triangle_is_line ())
    {
        std::pair<vector_t<TypeNum>, vector_t<TypeNum>> pair = select_ends_segment (
                                        other.get_a (), other.get_b (), other.get_c ());
        return check_intersection_tr_line (pair.first, pair.second);
    }

    if (!degenerate_tr () && other.triangle_is_point ())
    {
        return check_intersection_tr_point (other.get_a ());
    }

    return false;
}

template <typename TypeNum>
inline std::pair<vector_t<TypeNum>, vector_t<TypeNum>> triangle_t<TypeNum>::select_ends_segment (
             const vector_t<TypeNum>& p1, const vector_t<TypeNum>& p2, const vector_t<TypeNum>& p3)
{
    // c_ lies between a_ & b_
    if (check_intersection_line_point (p1, p2, p3)) 
        return {p1, p2};

    // b_ lies between a_ & c_
    else if (check_intersection_line_point (p1, p3, p2))
        return {p1, p3};

    // a_ lies between b_ & c_
    return {p2, p3};
};

template <typename TypeNum>
inline bool triangle_t<TypeNum>::check_intersection_tr_point (const vector_t<TypeNum>& p) const
{
    // guaranteed that point already lies in plane of triangle
    TypeNum res_1 = (vector_t{a_, b_}.cross_product(
                    vector_t{a_, p})).scalar_product(N_); 

    TypeNum res_2 = (vector_t{b_, c_}.cross_product(
                    vector_t{b_, p})).scalar_product(N_); 
    
    TypeNum res_3 = (vector_t{c_, a_}.cross_product(
                    vector_t{c_, p})).scalar_product(N_); 

    return ((res_1 >= -EPSILON<TypeNum> && res_2 >= -EPSILON<TypeNum> && res_3 >= -EPSILON<TypeNum>) || 
            (res_1 <= EPSILON<TypeNum> && res_2 <= EPSILON<TypeNum> && res_3 <= EPSILON<TypeNum>));
}

template <typename TypeNum>
inline bool triangle_t<TypeNum>::check_intersection_tr_line (const vector_t<TypeNum>& p1, 
                                                             const vector_t<TypeNum>& p2) const
{
    // guaranteed that line intersects plane of triangle (or lies)

    // lies
    if (N_.scalar_product(vector_t{p1, p2}) == 0)
    {
         return (check_intersection_line_line (p1, p2, a_, b_) ||
                 check_intersection_line_line (p1, p2, b_, c_) ||
                 check_intersection_line_line (p1, p2, c_, a_) ||
                 check_intersection_tr_point (p1) || check_intersection_tr_point (p2));
    }
    
    // intersects
    TypeNum t = - N_.scalar_product(vector_t{a_, p1}) /
                 N_.scalar_product(vector_t{p1, p2});

    vector_t p = p1 + ((p2 - p1) * t);

    return check_intersection_tr_point (p);
}

template <typename TypeNum>
inline bool triangle_t<TypeNum>::check_intersection_line_line (const vector_t<TypeNum>& line1_p1, 
                                                               const vector_t<TypeNum>& line1_p2, 
                                                               const vector_t<TypeNum>& line2_p1, 
                                                               const vector_t<TypeNum>& line2_p2)
{
    vector_t u{ line1_p1, line1_p2 }; // guiding vector line 1
    vector_t v{ line2_p1, line2_p2 }; // guiding vector line 2
    vector_t w{ line1_p1, line2_p1 }; // vector between line 1 & line 2

    // parallel lines
    if ((u.cross_product(v)).zero_vector()) 
    {
        return (check_intersection_line_point(line1_p1, line1_p2, line2_p1) || 
                check_intersection_line_point(line1_p1, line1_p2, line2_p2) ||
                check_intersection_line_point(line2_p1, line2_p2, line1_p1) || 
                check_intersection_line_point(line2_p1, line2_p2, line1_p2));
    }

    TypeNum u_u = u.scalar_product(u);
    TypeNum u_v = u.scalar_product(v);
    TypeNum v_v = v.scalar_product(v);
    TypeNum u_w = u.scalar_product(w);
    TypeNum v_w = v.scalar_product(w);

    TypeNum denominator = u_u * v_v - u_v * u_v;
    TypeNum parameter_1 = (v_v * u_w - u_v * v_w) / denominator;
    TypeNum parameter_2 = (u_v * u_w - u_u * v_w) / denominator;
  
    vector_t line1_parameter1 = line1_p1 + (line1_p2 - line1_p1) * parameter_1;
    vector_t line2_parameter2 = line2_p1 + (line2_p2 - line2_p1) * parameter_2;


    return ((vector_t{line1_parameter1, line2_parameter2}.zero_vector()) &&
            (parameter_1 >= -EPSILON<TypeNum>) && (parameter_1 <= 1 + EPSILON<TypeNum>)    &&
            (parameter_2 >= -EPSILON<TypeNum>) && (parameter_2 <= 1 + EPSILON<TypeNum>));
}

template <typename TypeNum>
inline bool triangle_t<TypeNum>::check_intersection_line_point (const vector_t<TypeNum>& line_p1, 
                                    const vector_t<TypeNum>& line_p2, const vector_t<TypeNum>& p)
{
    vector_t vec1{ line_p1, line_p2 };
    vector_t vec2{ line_p1, p };

    if (!vec1.cross_product (vec2).zero_vector ())
    {
        return false;
    }

    return ((std::fmin (line_p1.cor_x, line_p2.cor_x) - EPSILON<TypeNum> <= p.cor_x) &&
            (std::fmax (line_p1.cor_x, line_p2.cor_x) + EPSILON<TypeNum> >= p.cor_x) &&

            (std::fmin (line_p1.cor_y, line_p2.cor_y) - EPSILON<TypeNum> <= p.cor_y) &&
            (std::fmax (line_p1.cor_y, line_p2.cor_y) + EPSILON<TypeNum> >= p.cor_y) &&

            (std::fmin (line_p1.cor_z, line_p2.cor_z) - EPSILON<TypeNum> <= p.cor_z) &&
            (std::fmax (line_p1.cor_z, line_p2.cor_z) + EPSILON<TypeNum> >= p.cor_z));
}

template <typename TypeNum>
inline bool triangle_t<TypeNum>::check_intersection_point_point (const vector_t<TypeNum>& p1, 
                                                                 const vector_t<TypeNum>& p2)
{
    return (p1 == p2);
}

template <typename TypeNum>
inline bool triangle_t<TypeNum>::check_same_sign_distance (const triangle_t& other) const
{
    if (other.degenerate_tr())
    {
        return false;
    }
    TypeNum distance_1 = other.distance_point_plane_tr (a_);
    TypeNum distance_2 = other.distance_point_plane_tr (b_);
    TypeNum distance_3 = other.distance_point_plane_tr (c_);

    return ((distance_1 > EPSILON<TypeNum> && distance_2 > EPSILON<TypeNum> && distance_3 > EPSILON<TypeNum>) || 
            (distance_1 < -EPSILON<TypeNum> && distance_2 < -EPSILON<TypeNum> && distance_3 < -EPSILON<TypeNum>));
}

template <typename TypeNum>
inline bool triangle_t<TypeNum>::check_intersection_tr_of_line (const triangle_t& other) const
{
    // D = direction of the common line
    vector_t D = N_.cross_product (other.get_N ());

    // triangles lie in same plane
    if (D.zero_vector ())
    {
        return (check_intersection_tr_line (other.get_a (), other.get_b ()) ||
                check_intersection_tr_line (other.get_a (), other.get_c ()) ||
                check_intersection_tr_line (other.get_c (), other.get_b ()) ||

                other.check_intersection_tr_line (a_, b_) ||
                other.check_intersection_tr_line (a_, c_) ||
                other.check_intersection_tr_line (c_, b_));
    }

    TypeNum D_x = std::fabs (D.cor_x);
    TypeNum D_y = std::fabs (D.cor_y);
    TypeNum D_z = std::fabs (D.cor_z);

    axis_t axis = axis_t::Z_axis;

    if ((D_y - D_x) < EPSILON<TypeNum> && (D_z - D_x) < EPSILON<TypeNum>)
        axis = axis_t::X_axis;
    else if ((D_x - D_y) < EPSILON<TypeNum> && (D_z - D_y) < EPSILON<TypeNum>)
        axis = axis_t::Y_axis;

    std::pair<TypeNum, TypeNum> pair_1 = projection (axis, other);
    TypeNum t1 = pair_1.first;
    TypeNum t2 = pair_1.second;

    std::pair<TypeNum, TypeNum> pair_2 = other.projection (axis, *this);
    TypeNum t3 = pair_2.first;
    TypeNum t4 = pair_2.second;

    return ((std::min (t1, t2) <= std::max (t3, t4)) && 
            (std::min (t3, t4) <= std::max (t1, t2)));
}

template <typename TypeNum>
inline std::pair<TypeNum, TypeNum> triangle_t<TypeNum>::projection (axis_t axis, const triangle_t& other) const
{
    // point_1 and point_2 lie on the same side of tr_2 and that point_mid lies on
    // the other side
    vector_t point_1   = a_;
    vector_t point_mid = b_;
    vector_t point_2   = c_;

    // choose which axis to project on
    TypeNum project_point_1   = point_1.cor_z;
    TypeNum project_point_mid = point_mid.cor_z;
    TypeNum project_point_2   = point_2.cor_z;

    if (axis == axis_t::X_axis)
    {
        project_point_1   = point_1.cor_x;
        project_point_mid = point_mid.cor_x;
        project_point_2   = point_2.cor_x;
    }
    else if (axis == axis_t::Y_axis)
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
    TypeNum t1 = project_point_1 + 
                (project_point_mid - project_point_1) * (other.distance_point_plane_tr (point_1) / 
                (other.distance_point_plane_tr (point_1) - other.distance_point_plane_tr (point_mid)));

    TypeNum t2 = project_point_2 + 
                (project_point_mid - project_point_2) * (other.distance_point_plane_tr (point_2) / 
                (other.distance_point_plane_tr (point_2) - other.distance_point_plane_tr (point_mid)));

    return { t1, t2 };
}
} // ClassTriangle

// ----------------------------------------------------------------------------------

// ------------------------------OTHER_FUNC------------------------------------------

namespace Triangle {

template <typename TypeNum>
std::istream& operator>> (std::istream& in, vector_t<TypeNum>& p)
{
    in >> p.cor_x >> p.cor_y >> p.cor_z;
    return in;
}

inline void give_error_input_data ()
{
    std::cerr << "Invalid input" << std::endl;
    std::exit(EXIT_FAILURE);
}

} // Triangle

// ----------------------------------------------------------------------------------

#endif // TRIANGLES_HPP
