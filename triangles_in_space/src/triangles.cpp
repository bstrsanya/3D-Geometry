#include "triangles.hpp"

bool check_intersection (const triangle_t& tr_1, const triangle_t& tr_2)
{
    if (check_same_sign_distance (tr_1, tr_2) || 
        check_same_sign_distance (tr_2, tr_1))
    {
        return false; // one of the triangles lies completely in the half-plane of
                      // the other
    }

    // non-degenerate triangles
    if (!tr_1.degenerate_tr () && !tr_2.degenerate_tr ())
        return check_intersection_tr_of_line (tr_1, tr_2);

    // the same kind of degeneracy
    if (tr_1.triangle_is_line () && tr_2.triangle_is_line ())
    {
        std::pair<point_t, point_t> pair1 = select_ends_segment (
                                        tr_1.get_a (), tr_1.get_b (), tr_1.get_c ());
        std::pair<point_t, point_t> pair2 = select_ends_segment (
                                        tr_2.get_a (), tr_2.get_b (), tr_2.get_c ());
        return check_line_line (pair1.first, pair1.second, pair2.first, pair2.second);
    }

    if (tr_1.triangle_is_point () && tr_2.triangle_is_point ())
        return check_point_point (tr_1.get_a (), tr_2.get_a ());

    // different kinds of degeneracy
    return (check_different_degeneracies (tr_1, tr_2) || 
            check_different_degeneracies (tr_2, tr_1));
}

bool check_different_degeneracies (const triangle_t& tr_1, const triangle_t& tr_2)
{
    if (tr_1.triangle_is_line () && tr_2.triangle_is_point ())
    {
        std::pair<point_t, point_t> pair = select_ends_segment (
                                        tr_1.get_a (), tr_1.get_b (), tr_1.get_c ());
        return check_line_point (pair.first, pair.second, tr_2.get_a ());
    }

    if (!tr_1.degenerate_tr () && tr_2.triangle_is_line ())
    {
        std::pair<point_t, point_t> pair = select_ends_segment (
                                        tr_2.get_a (), tr_2.get_b (), tr_2.get_c ());
        return check_triangle_line (tr_1, pair.first, pair.second);
    }

    if (!tr_1.degenerate_tr () && tr_2.triangle_is_point ())
    {
        return check_triangle_point (tr_1, tr_2.get_a ());
    }

    return false;
}

std::pair<point_t, point_t> select_ends_segment (const point_t& p1, 
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

bool check_triangle_point (const triangle_t& tr_1, const point_t& p)
{
    // guaranteed that point already lies in plane of triangle
    double res_1 = (vector_t{tr_1.get_a(), tr_1.get_b()}.cross_product(
                    vector_t{tr_1.get_a(), p})).scalar_product(tr_1.get_N()); 

    double res_2 = (vector_t{tr_1.get_b(), tr_1.get_c()}.cross_product(
                    vector_t{tr_1.get_b(), p})).scalar_product(tr_1.get_N()); 
    
    double res_3 = (vector_t{tr_1.get_c(), tr_1.get_a()}.cross_product(
                    vector_t{tr_1.get_c(), p})).scalar_product(tr_1.get_N()); 

    return ((res_1 >= -epsilon && res_2 >= -epsilon && res_3 >= -epsilon) || 
            (res_1 <= epsilon && res_2 <= epsilon && res_3 <= epsilon));
}

bool check_triangle_line (const triangle_t& tr_1, const point_t& p1, const point_t& p2)
{
    // guaranteed that line intersects plane of triangle (or lies)

    // lies
    if (tr_1.get_N().scalar_product(vector_t{p1, p2}) == 0)
    {
         return (check_line_line (p1, p2, tr_1.get_a(), tr_1.get_b()) ||
                 check_line_line (p1, p2, tr_1.get_b(), tr_1.get_c()) ||
                 check_line_line (p1, p2, tr_1.get_c(), tr_1.get_a()) ||
                 check_triangle_point (tr_1, p1) || check_triangle_point (tr_1, p2));
    }
    
    // intersects
    double t = - tr_1.get_N().scalar_product(vector_t{tr_1.get_a(), p1}) /
                 tr_1.get_N().scalar_product(vector_t{p1, p2});

    point_t p = p1 + ((p2 - p1) * t);

    return check_triangle_point (tr_1, p);
}

bool check_line_line (const point_t& line1_p1, const point_t& line1_p2,
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
            (parameter_1 >= -epsilon) && (parameter_1 <= 1 + epsilon)    &&
            (parameter_2 >= -epsilon) && (parameter_2 <= 1 + epsilon));
}

bool check_line_point (const point_t& line_p1, const point_t& line_p2, const point_t& p)
{
    vector_t vec1{ line_p1, line_p2 };
    vector_t vec2{ line_p1, p };

    if (!vec1.cross_product (vec2).zero_vector ())
    {
        return false;
    }

    return ((std::fmin (line_p1.get_x (), line_p2.get_x ()) - epsilon <= p.get_x ()) &&
            (std::fmax (line_p1.get_x (), line_p2.get_x ()) + epsilon >= p.get_x ()) &&

            (std::fmin (line_p1.get_y (), line_p2.get_y ()) - epsilon <= p.get_y ()) &&
            (std::fmax (line_p1.get_y (), line_p2.get_y ()) + epsilon >= p.get_y ()) &&

            (std::fmin (line_p1.get_z (), line_p2.get_z ()) - epsilon <= p.get_z ()) &&
            (std::fmax (line_p1.get_z (), line_p2.get_z ()) + epsilon >= p.get_z ()));
}

bool check_point_point (const point_t& p1, const point_t& p2)
{
    return (p1 == p2);
}

bool check_same_sign_distance (const triangle_t& tr_1, const triangle_t& tr_2)
{
    if (tr_2.degenerate_tr())
    {
        return false;
    }
    double distance_1 = tr_2.distance_point_plane_tr (tr_1.get_a ());
    double distance_2 = tr_2.distance_point_plane_tr (tr_1.get_b ());
    double distance_3 = tr_2.distance_point_plane_tr (tr_1.get_c ());

    return ((distance_1 > epsilon && distance_2 > epsilon && distance_3 > epsilon) || 
            (distance_1 < -epsilon && distance_2 < -epsilon && distance_3 < -epsilon));
}

bool check_intersection_tr_of_line (const triangle_t& tr_1, const triangle_t& tr_2)
{
    // D = direction of the common line
    vector_t D = tr_1.get_N ().cross_product (tr_2.get_N ());

    // triangles lie in same plane
    if (D.zero_vector ())
    {
        return (check_triangle_line (tr_1, tr_2.get_a (), tr_2.get_b ()) ||
                check_triangle_line (tr_1, tr_2.get_a (), tr_2.get_c ()) ||
                check_triangle_line (tr_1, tr_2.get_c (), tr_2.get_b ()) ||

                check_triangle_point (tr_1, tr_2.get_a ()) ||
                check_triangle_point (tr_1, tr_2.get_b ()) || 
                check_triangle_point (tr_1, tr_2.get_c ()) ||
                
                check_triangle_line (tr_2, tr_1.get_a (), tr_1.get_b ()) ||
                check_triangle_line (tr_2, tr_1.get_a (), tr_1.get_c ()) ||
                check_triangle_line (tr_2, tr_1.get_c (), tr_1.get_b ()) ||

                check_triangle_point (tr_2, tr_1.get_a ()) ||
                check_triangle_point (tr_2, tr_1.get_b ()) || 
                check_triangle_point (tr_2, tr_1.get_c ()));
    }

    double D_x = std::fabs (D.get_x ());
    double D_y = std::fabs (D.get_y ());
    double D_z = std::fabs (D.get_z ());

    char axis = 'z';

    if ((D_y - D_x) < epsilon && (D_z - D_x) < epsilon)
        axis = 'x';
    else if ((D_x - D_y) < epsilon && (D_z - D_y) < epsilon)
        axis = 'y';

    std::pair<double, double> pair_1 = projection (axis, tr_1, tr_2);
    double t1 = pair_1.first;
    double t2 = pair_1.second;

    std::pair<double, double> pair_2 = projection (axis, tr_2, tr_1);
    double t3 = pair_2.first;
    double t4 = pair_2.second;

    return ((std::min (t1, t2) <= std::max (t3, t4)) && 
            (std::min (t3, t4) <= std::max (t1, t2)));
}

std::pair<double, double> projection (char axis, const triangle_t& tr_1, const triangle_t& tr_2)
{
    // point_1 and point_2 lie on the same side of tr_2 and that point_mid lies on
    // the other side
    point_t point_1   = tr_1.get_a ();
    point_t point_mid = tr_1.get_b ();
    point_t point_2   = tr_1.get_c ();

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
    if (tr_2.point_lie_in_plane_tr (tr_1.get_a ()) && 
        tr_2.point_lie_in_plane_tr (tr_1.get_b ()))
    {
        return { project_point_1, project_point_mid };
    }

    if (tr_2.point_lie_in_plane_tr (tr_1.get_b ()) && 
        tr_2.point_lie_in_plane_tr (tr_1.get_c ()))
    {
        return { project_point_2, project_point_mid };
    }
    if (tr_2.point_lie_in_plane_tr (tr_1.get_a ()) && 
        tr_2.point_lie_in_plane_tr (tr_1.get_c ()))
    {
        return { project_point_1, project_point_2 };
    }

    // By default, the middle point is point B, below are 3 cases where point B cannot be middle

    if (tr_2.distance_point_plane_tr (tr_1.get_b ()) * 
        tr_2.distance_point_plane_tr (tr_1.get_c ()) > 0)
    {
        std::swap (project_point_1, project_point_mid);
        std::swap (point_1, point_mid);
    }
    else if (tr_2.distance_point_plane_tr (tr_1.get_a ()) * 
             tr_2.distance_point_plane_tr (tr_1.get_b ()) > 0)
    {
        std::swap (project_point_2, project_point_mid);
        std::swap (point_2, point_mid);
    }

    if (tr_2.point_lie_in_plane_tr (tr_1.get_b ()) && 
        tr_2.distance_point_plane_tr (tr_1.get_a ()) * 
        tr_2.distance_point_plane_tr (tr_1.get_c ()) < 0)
    {
        std::swap (project_point_1, project_point_mid);
        std::swap (point_1, point_mid);
    }

    // counting projections on the selected axis
    double t1 = project_point_1 + 
                (project_point_mid - project_point_1) * (tr_2.distance_point_plane_tr (point_1) / 
                (tr_2.distance_point_plane_tr (point_1) - tr_2.distance_point_plane_tr (point_mid)));

    double t2 = project_point_2 + 
                (project_point_mid - project_point_2) * (tr_2.distance_point_plane_tr (point_2) / 
                (tr_2.distance_point_plane_tr (point_2) - tr_2.distance_point_plane_tr (point_mid)));

    return { t1, t2 };
}

std::istream& operator>> (std::istream& in, point_t& p)
{
    in >> p.x_ >> p.y_ >> p.z_;
    return in;
}