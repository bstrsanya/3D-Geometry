#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

#include "triangles.hpp"
#include "grid.hpp"

int main ()
{
    std::size_t N = 0;
    if (!(std::cin >> N) || N > 1'000'000)
    {
        std::cerr << "The number of triangles should be in [0, 1'000'000]" << std::endl;
        return 1;
    }

    if (N == 0) 
    {
        return 0;
    }

    std::vector<Triangle::triangle_t<float>> array_triangle(N);

    for (std::size_t i = 0; i < N; ++i)
    {
        Triangle::vector_t<float> p1, p2, p3;
        if (!(std::cin >> p1 >> p2 >> p3)) 
        {
            std::cerr << "Expected (<number_of_triangles> * 9) numerical values" << std::endl;
            return 1;
        }
        array_triangle[i] = { p1, p2, p3 };
    }

    Triangle::adaptive_grid_t tree(std::move (array_triangle));
    auto triangle_num = tree.get_num_tr_intersection ();
    
    for (auto tmp: triangle_num)
    {
         std::cout << tmp << "\n";
    }
}
