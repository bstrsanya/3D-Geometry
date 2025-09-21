#include "triangles.hpp"
#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

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

    std::set<std::size_t> triangle_num;

    for (std::size_t i = 0; i < N; ++i)
    {
        for (std::size_t j = i + 1; j < N; ++j)
        {
            if (check_intersection (array_triangle[i], array_triangle[j]))
            {
                triangle_num.insert (i);
                triangle_num.insert (j);
            }
        }
    }

    for (auto tmp: triangle_num)
    {
        std::cout << tmp + 1 << "\n";
    }
}
