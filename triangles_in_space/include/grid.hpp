#ifndef ADAPTIVE_GRID_HPP
#define ADAPTIVE_GRID_HPP

#include <set>
#include <vector>
#include <array>
#include <limits>

#include "triangles.hpp"

// --------------------------------GRID_T--------------------------------------------
namespace Triangle {

template <typename TypeNum>
class adaptive_grid_t
{
private:
    static const std::size_t OPTIMAL_NUM_TR_IN_SPACE = 15;
    static const std::size_t MAX_VALUE_DEEP_RECURSION = 6;
    static const std::size_t OCTREE_CHILD_COUNT = 8;

    std::vector<triangle_t<TypeNum>> array_triangle_;
    std::vector<std::vector<std::size_t>> array_leaf_tree_ {};
    std::set<std::size_t> num_tr_intersection_ {};   

    TypeNum count_bounding_cube ();
    TypeNum nearest_power_of_two (TypeNum num);
    void recursive_construction_grid (const vector_t<TypeNum>& p_min, 
                                      const vector_t<TypeNum>& p_max,
                                      std::vector<std::size_t>& num_triangles, int dep);
    void naive_verification (std::vector<std::size_t>& num1);  

public:
    adaptive_grid_t (std::vector<triangle_t<TypeNum>>&& array_triangle);
    const std::set<std::size_t>& get_num_tr_intersection ();
};

template <typename TypeNum>
inline adaptive_grid_t<TypeNum>::adaptive_grid_t (std::vector<triangle_t<TypeNum>>&& array_triangle)
                                                  : array_triangle_ (std::move (array_triangle))
{
    TypeNum max_coordinate = count_bounding_cube ();

    vector_t p_max = vector_t (max_coordinate, max_coordinate, max_coordinate);
    vector_t p_min = vector_t (-max_coordinate, -max_coordinate, -max_coordinate);

    std::vector<std::size_t> num_triangles{};
    for (std::size_t i = 0; i < array_triangle_.size(); i++)
    {
        num_triangles.push_back (i);
    }

    recursive_construction_grid (p_min, p_max, num_triangles, 0);
}

template <typename TypeNum>
inline TypeNum adaptive_grid_t<TypeNum>::count_bounding_cube ()
{
    // definition of a symmetric cube (relative to the point O(0,0,0)) that contains all triangles
    vector_t<TypeNum> p_min = array_triangle_[0].get_a();
    vector_t<TypeNum> p_max = array_triangle_[0].get_a();
    for (const auto &tr : array_triangle_) 
    {
        const vector_t<TypeNum> verts[] = { tr.get_a(), tr.get_b(), tr.get_c() };
        for (const auto &v : verts)
        {
            p_min.cor_x = std::min(p_min.cor_x, v.cor_x);
            p_min.cor_y = std::min(p_min.cor_y, v.cor_y);
            p_min.cor_z = std::min(p_min.cor_z, v.cor_z);

            p_max.cor_y = std::max(p_max.cor_y, v.cor_y);
            p_max.cor_x = std::max(p_max.cor_x, v.cor_x);
            p_max.cor_z = std::max(p_max.cor_z, v.cor_z);
        }
    }

    TypeNum max_coordinate = std::max (std::fabs (p_min.cor_x), std::fabs (p_max.cor_x));
    max_coordinate = std::max (max_coordinate, 
                               std::max (std::fabs (p_min.cor_y), std::fabs (p_max.cor_y)));
    max_coordinate = std::max (max_coordinate, 
                               std::max (std::fabs (p_min.cor_z), std::fabs (p_max.cor_z)));

    max_coordinate = nearest_power_of_two (max_coordinate);

    return max_coordinate;
}

template <typename TypeNum>
inline TypeNum adaptive_grid_t<TypeNum>::nearest_power_of_two (TypeNum num)
{
    if (num <= static_cast<TypeNum>(0))
    {
        return static_cast<TypeNum>(1);
    }

    TypeNum exp = std::ceil(std::log2(num));
    TypeNum result = std::pow(static_cast<TypeNum>(2), exp);
    return result;
}

template <typename TypeNum>
inline void adaptive_grid_t<TypeNum>::recursive_construction_grid (const vector_t<TypeNum>& p_min, 
                                                                   const vector_t<TypeNum>& p_max,
                                     std::vector<std::size_t>& num_triangles, int depth_recursion)
{
    depth_recursion++;
    
    // recursion exit conditions
    if (num_triangles.size () <= OPTIMAL_NUM_TR_IN_SPACE || 
               depth_recursion > MAX_VALUE_DEEP_RECURSION)
    {
        array_leaf_tree_.push_back (num_triangles);
        return;
    }

    // common point
    vector_t central_point = (p_min + p_max) / 2;

    // dividing space into 8 equal parts that have a common point
    std::array<std::vector<std::size_t>, OCTREE_CHILD_COUNT> array_space{};
    std::array<vector_t<TypeNum>, OCTREE_CHILD_COUNT> array_point = {
                                          vector_t {p_max.cor_x, p_max.cor_y, p_max.cor_z},
                                          vector_t {p_min.cor_x, p_max.cor_y, p_max.cor_z},
                                          vector_t {p_min.cor_x, p_min.cor_y, p_max.cor_z},
                                          vector_t {p_max.cor_x, p_min.cor_y, p_max.cor_z},
                                          vector_t {p_max.cor_x, p_max.cor_y, p_min.cor_z},
                                          vector_t {p_min.cor_x, p_max.cor_y, p_min.cor_z},
                                          vector_t {p_min.cor_x, p_min.cor_y, p_min.cor_z},
                                          vector_t {p_max.cor_x, p_min.cor_y, p_min.cor_z}};


    // distribution of triangles into new subspaces
    for (const auto& n_tr : num_triangles)
    {
        triangle_t<TypeNum>& tr = array_triangle_[n_tr];
        for (std::size_t i = 0; i < OCTREE_CHILD_COUNT; i++)
        {
            if (tr.triangle_lie_in_space (central_point, array_point[i]))
            {
                array_space[i].push_back (n_tr);
            }
        }
    }

    // for non-empty subspaces we call a recursive partition
    for (std::size_t i = 0; i < OCTREE_CHILD_COUNT; i++)
    {
        if (!array_space[i].empty())
        {
            recursive_construction_grid (central_point, array_point[i], 
                                         array_space[i], depth_recursion);
        }
    }
}

template <typename TypeNum>
inline const std::set<std::size_t>& adaptive_grid_t<TypeNum>::get_num_tr_intersection ()
{
    // naive intersection count for each part of the partition
    for (auto& leaf : array_leaf_tree_)
    {
        naive_verification (leaf);
    }
    return num_tr_intersection_;
}

template <typename TypeNum>
inline void adaptive_grid_t<TypeNum>::naive_verification (std::vector<std::size_t>& num)
{
    for (auto it1 = num.begin(); it1 != num.end(); ++it1)
    {
        auto it2 = it1;
        ++it2;

        for (; it2 != num.end(); ++it2)
        {
            if (array_triangle_[*it1].check_intersection_tr_tr(array_triangle_[*it2]))
            {
                num_tr_intersection_.insert (*it1);
                num_tr_intersection_.insert (*it2);
            }
        }
    }
}
} // ClassAdaptiveGrid

// ----------------------------------------------------------------------------------

#endif // ADAPTIVE_GRID_HPP