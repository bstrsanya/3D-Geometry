#include "triangles.hpp"

bool check_intersection (const triangle_t& tr_1, const triangle_t& tr_2)
{
    if (check_same_sign_distance (tr_1, tr_2) || 
        check_same_sign_distance (tr_2, tr_1))
    {
        return false; // one of the triangles lies completely in the half-plane of
                      // the other
    }

    return check_intersection_tr_of_line (tr_1, tr_2);
}

bool check_same_sign_distance (const triangle_t& tr_1, const triangle_t& tr_2)
{
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

    // otherwise, looking for the right point_mid
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