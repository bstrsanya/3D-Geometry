#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

#include "triangles.hpp"
#include "grid.hpp"

int main ()
{
    std::size_t N = 0;
    std::cin >> N;
    std::vector<triangle_t<float>> array_triangle(N);

    for (std::size_t i = 0; i < N; ++i)
    {
        vector_t<float> p1, p2, p3;
        if (!(std::cin >> p1 >> p2 >> p3)) 
        {
            std::cerr << "Invalid input" << std::endl;
            std::exit(EXIT_FAILURE);
        }
        array_triangle[i] = { p1, p2, p3 };
    }

    adaptive_grid_t tree(std::move (array_triangle));
    std::set<std::size_t> triangle_num = tree.get_num_tr_intersection ();
    
    for (auto tmp: triangle_num)
    {
         std::cout << tmp << "\n";
    }
}
