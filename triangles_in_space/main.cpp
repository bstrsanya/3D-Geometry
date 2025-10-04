#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

#include "triangles.hpp"
#include "octree.hpp"

int main ()
{
    std::size_t N = 0;
    std::cin >> N;
    std::vector<triangle_t> array_triangle;

    for (std::size_t i = 0; i < N; ++i)
    {
        point_t p1, p2, p3;
        std::cin >> p1 >> p2 >> p3;
        array_triangle.push_back ({ p1, p2, p3 });
    }

    octree_t tree(array_triangle);
    std::set<std::size_t> triangle_num = tree.get_num_tr_intersection ();
    
    for (auto tmp: triangle_num)
    {
         std::cout << tmp << "\n";
    }
}
