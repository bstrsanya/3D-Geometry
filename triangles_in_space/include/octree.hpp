#ifndef OCTREE_HPP
#define OCTREE_HPP

#include <set>
#include <vector>
#include <array>
#include <limits>

#include "triangles.hpp"

const double MAX_DOUBLE = std::numeric_limits<double>::max();
const double MIN_DOUBLE = -std::numeric_limits<double>::max();
const std::size_t OPTIMAL_NUM_TR_IN_SPACE = 15;
const std::size_t MAX_VALUE_DEEP_RECURSION = 6;
const std::size_t OCTREE_CHILD_COUNT = 8;

// ------------------------------NODE_T----------------------------------------------

class node_t
{
private:
    std::vector<std::size_t> num_triangles_in_same_space_ {};

public:
    node_t (const std::vector<std::size_t>& num_tr) : num_triangles_in_same_space_(num_tr) {};
    std::vector<std::size_t>& get_num_triangles () { return num_triangles_in_same_space_; }
};

// ----------------------------------------------------------------------------------

// ------------------------------OCTREE_T--------------------------------------------

class octree_t
{
private:
    std::vector<triangle_t>& array_triangle_;
    std::vector<node_t*> array_leaf_tree_ {};
    std::set<std::size_t> num_tr_intersection_ {};   

    double count_bounding_cube ();
    double nearest_power_of_two (double num);
    void recursive_construction_tree (const point_t& p_min, const point_t& p_max,
                                      std::vector<std::size_t>& num_triangles, int dep);

    void naive_verification (std::vector<std::size_t>& num1);  

public:
    octree_t (std::vector<triangle_t>& array_triangle);
    ~octree_t() { for (auto& tmp : array_leaf_tree_) { delete tmp; } };

    std::set<std::size_t> get_num_tr_intersection ();
};

inline octree_t::octree_t (std::vector<triangle_t>& array_triangle) : array_triangle_(array_triangle)
{
    double max_coordinate = count_bounding_cube ();

    point_t p_max = point_t (max_coordinate, max_coordinate, max_coordinate);
    point_t p_min = point_t (-max_coordinate, -max_coordinate, -max_coordinate);

    std::vector<std::size_t> num_triangles{};
    for (std::size_t i = 0; i < array_triangle.size(); i++)
    {
        num_triangles.push_back (i);
    }

    recursive_construction_tree (p_min, p_max, num_triangles, 0);
}

inline double octree_t::count_bounding_cube ()
{
    point_t p_min {MAX_DOUBLE, MAX_DOUBLE, MAX_DOUBLE};
    point_t p_max {MIN_DOUBLE, MIN_DOUBLE, MIN_DOUBLE};
    for (const auto &tr : array_triangle_) 
    {
        const point_t verts[] = { tr.get_a(), tr.get_b(), tr.get_c() };
        for (const auto &v : verts)
        {
            p_min.x_ = std::min(p_min.x_, v.x_);
            p_min.y_ = std::min(p_min.y_, v.y_);
            p_min.z_ = std::min(p_min.z_, v.z_);

            p_max.y_ = std::max(p_max.y_, v.y_);
            p_max.x_ = std::max(p_max.x_, v.x_);
            p_max.z_ = std::max(p_max.z_, v.z_);
        }
    }

    double max_coordinate = std::max (std::fabs (p_min.x_), std::fabs (p_max.x_));
    max_coordinate = std::max (max_coordinate, std::max (std::fabs (p_min.y_), std::fabs (p_max.y_)));
    max_coordinate = std::max (max_coordinate, std::max (std::fabs (p_min.z_), std::fabs (p_max.z_)));

    max_coordinate = nearest_power_of_two (max_coordinate);

    return max_coordinate;
}

inline double octree_t::nearest_power_of_two (double num)
{
    int x = static_cast<int> (num) + 1; 
    if (x == 0) 
        return static_cast<double> (1);
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return static_cast<double> (x + 1);
}

inline void octree_t::recursive_construction_tree (const point_t& p_min, const point_t& p_max,
                                         std::vector<std::size_t>& num_triangles, int depth_recursion)
{
    depth_recursion++;
    
    if (num_triangles.size () <= OPTIMAL_NUM_TR_IN_SPACE || 
               depth_recursion > MAX_VALUE_DEEP_RECURSION)
    {
        node_t* main_node = new node_t{num_triangles};
        array_leaf_tree_.push_back (main_node);
        return;
    }

    point_t central_point = (p_min + p_max) / 2;

    std::array<std::vector<std::size_t>, OCTREE_CHILD_COUNT> array_space{};
    std::array<point_t, OCTREE_CHILD_COUNT> array_point = {
                                          point_t {p_max.x_, p_max.y_, p_max.z_},
                                          point_t {p_min.x_, p_max.y_, p_max.z_},
                                          point_t {p_min.x_, p_min.y_, p_max.z_},
                                          point_t {p_max.x_, p_min.y_, p_max.z_},
                                          point_t {p_max.x_, p_max.y_, p_min.z_},
                                          point_t {p_min.x_, p_max.y_, p_min.z_},
                                          point_t {p_min.x_, p_min.y_, p_min.z_},
                                          point_t {p_max.x_, p_min.y_, p_min.z_}};                                          ;


    for (auto n_tr = num_triangles.begin(); n_tr != num_triangles.end(); n_tr++)
    {
        triangle_t& tr = array_triangle_[*n_tr];
        for (std::size_t i = 0; i < OCTREE_CHILD_COUNT; i++)
        {
            if (tr.triangle_lie_in_space (central_point, array_point[i]))
            {
                array_space[i].push_back (*n_tr);
            }
        }
    }

    for (std::size_t i = 0; i < OCTREE_CHILD_COUNT; i++)
    {
        if (!array_space[i].empty())
        {
            recursive_construction_tree (central_point, array_point[i], array_space[i], depth_recursion);
        }
    }
}

inline std::set<std::size_t> octree_t::get_num_tr_intersection ()
{
    for (auto& leaf : array_leaf_tree_)
    {
        naive_verification (leaf->get_num_triangles ());
    }
    return num_tr_intersection_;
}

inline void octree_t::naive_verification (std::vector<std::size_t>& num)
{
    for (auto it1 = num.begin(); it1 != num.end(); ++it1)
    {
        auto it2 = it1;
        ++it2;

        for (; it2 != num.end(); ++it2)
        {
            if (array_triangle_[*it1].check_intersection(array_triangle_[*it2]))
            {
                num_tr_intersection_.insert (*it1);
                num_tr_intersection_.insert (*it2);
            }
        }
    }
}

// ----------------------------------------------------------------------------------

#endif // OCTREE_HPP