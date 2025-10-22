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
        Triangle::give_error_input_data ();
    }
    std::vector<Triangle::triangle_t<float>> array_triangle(N);

    for (std::size_t i = 0; i < N; ++i)
    {
        Triangle::vector_t<float> p1, p2, p3;
        if (!(std::cin >> p1 >> p2 >> p3)) 
        {
            Triangle::give_error_input_data ();
        }
        array_triangle[i] = { p1, p2, p3 };
    }

    Triangle::adaptive_grid_t tree(std::move (array_triangle));
    std::set<std::size_t> triangle_num = tree.get_num_tr_intersection ();
    
    for (auto tmp: triangle_num)
    {
         std::cout << tmp << "\n";
    }
}
