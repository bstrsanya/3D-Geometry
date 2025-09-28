#include <gtest/gtest.h>

#include "./../include/triangles.hpp"

// ------------------------------TESTING_SCALAR_PRODUCT------------------------------

TEST (test_vector, scalar_product_easy)
{
    vector_t vec_1 { 1.0, 1.5, 2.0 };
    vector_t vec_2 { 3.0, 5.0, 4.4 };
    double answer = vec_1.scalar_product (vec_2);
    EXPECT_DOUBLE_EQ (answer, 19.3);
}

TEST (test_vector, scalar_product_zero)
{
    vector_t vec_1 { 0.0, 0.0, 0.0 };
    vector_t vec_2 { 3.0, -5.0, 10.0 };
    double answer = vec_1.scalar_product (vec_2);
    EXPECT_DOUBLE_EQ (answer, 0.0);
}

TEST (test_vector, scalar_product_orthogonal)
{
    vector_t vec_1 { 1.0, 0.0, 0.0 };
    vector_t vec_2 { 0.0, 1.0, 0.0 };
    double answer = vec_1.scalar_product (vec_2);
    EXPECT_DOUBLE_EQ (answer, 0.0);
}

TEST (test_vector, scalar_product_parallel)
{
    vector_t vec_1 { 2.0, 2.0, 2.0 };
    vector_t vec_2 { 1.0, 1.0, 1.0 };
    double answer = vec_1.scalar_product (vec_2);
    EXPECT_DOUBLE_EQ (answer, 6.0);
}

TEST (test_vector, scalar_product_large)
{
    vector_t vec_1 { 1e6, 2e6, -3e6 };
    vector_t vec_2 { -4e6, 5e6, 6e6 };
    double answer = vec_1.scalar_product (vec_2);
    EXPECT_DOUBLE_EQ (answer, -1.2e13);
}

TEST (test_vector, scalar_product_mixed)
{
    vector_t vec_1 { 0.1, 0.01, -0.001 };
    vector_t vec_2 { 1000.0, -100.0, 10.0 };
    double answer = vec_1.scalar_product (vec_2);
    EXPECT_DOUBLE_EQ (answer, 98.99);
}

TEST (test_vector, scalar_product_negative)
{
    vector_t vec_1 { -1.1, -2.2, -3.3 };
    vector_t vec_2 { -4.4, -5.5, -6.6 };
    double answer = vec_1.scalar_product (vec_2);
    EXPECT_DOUBLE_EQ (answer, 38.72);
}

// ----------------------------------------------------------------------------------

// ------------------------------TESTING_CROSS_PRODUCT-------------------------------

TEST (test_vector, cross_product_easy)
{
    vector_t vec_1 { 1.0, 1.5, 2.0 };
    vector_t vec_2 { 3.0, 5.0, 4.4 };
    vector_t answer = vec_1.cross_product (vec_2);
    EXPECT_DOUBLE_EQ (answer.get_x (), -3.4);
    EXPECT_DOUBLE_EQ (answer.get_y (), 1.6);
    EXPECT_DOUBLE_EQ (answer.get_z (), 0.5);
}

TEST (test_vector, cross_product_zero)
{
    vector_t vec_1 { 0.0, 0.0, 0.0 };
    vector_t vec_2 { 3.0, -5.0, 10.0 };
    vector_t answer = vec_1.cross_product (vec_2);
    EXPECT_DOUBLE_EQ (answer.get_x (), 0.0);
    EXPECT_DOUBLE_EQ (answer.get_y (), 0.0);
    EXPECT_DOUBLE_EQ (answer.get_z (), 0.0);
}

TEST (test_vector, cross_product_parallel)
{
    vector_t vec_1 { 2.0, 2.0, 2.0 };
    vector_t vec_2 { 4.0, 4.0, 4.0 };
    vector_t answer = vec_1.cross_product (vec_2);
    EXPECT_DOUBLE_EQ (answer.get_x (), 0.0);
    EXPECT_DOUBLE_EQ (answer.get_y (), 0.0);
    EXPECT_DOUBLE_EQ (answer.get_z (), 0.0);
}

TEST (test_vector, cross_product_orthogonal)
{
    vector_t vec_1 { 1.0, 0.0, 0.0 };
    vector_t vec_2 { 0.0, 1.0, 0.0 };
    vector_t answer = vec_1.cross_product (vec_2);
    EXPECT_DOUBLE_EQ (answer.get_x (), 0.0);
    EXPECT_DOUBLE_EQ (answer.get_y (), 0.0);
    EXPECT_DOUBLE_EQ (answer.get_z (), 1.0);
}

TEST (test_vector, cross_product_middle)
{
    vector_t vec_1 { 1.25, -3.5, 4.75 };
    vector_t vec_2 { -2.0, 0.5, 3.2 };
    vector_t answer = vec_1.cross_product (vec_2);
    EXPECT_DOUBLE_EQ (answer.get_x (), -13.575);
    EXPECT_DOUBLE_EQ (answer.get_y (), -13.5);
    EXPECT_DOUBLE_EQ (answer.get_z (), -6.375);
}

TEST (test_vector, cross_product_large)
{
    vector_t vec_1 { 1e6, 2e6, -3e6 };
    vector_t vec_2 { -4e6, 5e6, 6e6 };
    vector_t answer = vec_1.cross_product (vec_2);
    EXPECT_DOUBLE_EQ (answer.get_x (), 2.7e13);
    EXPECT_DOUBLE_EQ (answer.get_y (), 6e12);
    EXPECT_DOUBLE_EQ (answer.get_z (), 1.3e13);
}

TEST (test_vector, cross_product_small)
{
    vector_t vec_1 { 0.001, 0.002, 0.003 };
    vector_t vec_2 { 0.004, 0.005, 0.006 };
    vector_t answer = vec_1.cross_product (vec_2);
    EXPECT_DOUBLE_EQ (answer.get_x (), -0.000003);
    EXPECT_DOUBLE_EQ (answer.get_y (), 0.000006);
    EXPECT_DOUBLE_EQ (answer.get_z (), -0.000003);
}

TEST (test_vector, cross_product_negative)
{
    vector_t vec_1 { -1.1, -2.2, -3.3 };
    vector_t vec_2 { -4.4, -5.5, -6.6 };
    vector_t answer = vec_1.cross_product (vec_2);
    EXPECT_DOUBLE_EQ (answer.get_x (), -3.63);
    EXPECT_DOUBLE_EQ (answer.get_y (), 7.26);
    EXPECT_DOUBLE_EQ (answer.get_z (), -3.63);
}

// ----------------------------------------------------------------------------------

// ------------------------------TESTING_POINT_LIE_IN_PLANE_TR-----------------------

TEST (test_triangle, easy_1)
{
    point_t t1 (5.0, 0.0, 2.7);
    point_t t2 (1.0, 0.0, 10.3);
    point_t t3 (2.0, 0.0, 1.0);
    triangle_t tr (t1, t2, t3);
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 13.0, 0.0, 25.9 }));
    EXPECT_FALSE (tr.point_lie_in_plane_tr ({ 1.0, 1.0, 1.0 }));
}

TEST (test_triangle, easy_2)
{
    point_t t1 (8.9, -7.4, -10.25);
    point_t t2 (10.0, 0.0, -3.8);
    point_t t3 (-2.0, 16.0, 17.16);
    triangle_t tr (t1, t2, t3);
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 8.9, -7.4, -10.25 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 10.0, 0.0, -3.8 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ -2.0, 16.0, 17.16 }));
    EXPECT_FALSE (tr.point_lie_in_plane_tr ({ 1.0, 1.0, 1.0 }));
}

TEST (test_triangle, easy_3)
{
    point_t t1 (0.0, 0.0, 0.0);
    point_t t2 (1.0, 0.0, 0.0);
    point_t t3 (0.0, 1.0, 0.0);
    triangle_t tr (t1, t2, t3);

    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 0.0, 0.0, 0.0 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 1.0, 0.0, 0.0 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 0.0, 1.0, 0.0 }));

    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 0.5, 0.5, 0.0 }));

    EXPECT_FALSE (tr.point_lie_in_plane_tr ({ 0.5, 0.5, 1.0 }));
}

TEST (test_triangle, hard_1)
{
    point_t t1 (1.0, 2.0, 3.0);
    point_t t2 (-4.0, 5.5, 6.6);
    point_t t3 (7.7, -8.8, 9.9);
    triangle_t tr (t1, t2, t3);

    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 1.0, 2.0, 3.0 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ -4.0, 5.5, 6.6 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 7.7, -8.8, 9.9 }));

    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ (1.0 - 4.0 + 7.7) / 3.0,
                                             (2.0 + 5.5 - 8.8) / 3.0,
                                             (3.0 + 6.6 + 9.9) / 3.0 }));

    EXPECT_FALSE (tr.point_lie_in_plane_tr ({ 10.0, 10.0, 10.0 }));
}

TEST (test_triangle, hard_2)
{
    point_t t1 (-2.5, 4.1, 0.0);
    point_t t2 (3.3, -1.7, 0.0);
    point_t t3 (0.0, 0.0, 0.0);
    triangle_t tr (t1, t2, t3);

    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ -2.5, 4.1, 0.0 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 3.3, -1.7, 0.0 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 0.0, 0.0, 0.0 }));

    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 100.0, -200.0, 0.0 }));

    EXPECT_FALSE (tr.point_lie_in_plane_tr ({ 0.0, 0.0, 1.0 }));
}

TEST (test_triangle, hard_large)
{
    point_t t1 (1e8, 2e8, -3e8);
    point_t t2 (-4e8, 5e8, 6e8);
    point_t t3 (7e8, -8e8, 9e8);
    triangle_t tr (t1, t2, t3);

    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 1e8, 2e8, -3e8 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ -4e8, 5e8, 6e8 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 7e8, -8e8, 9e8 }));

    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ (1e8 - 4e8 + 7e8) / 3.0,
                                             (2e8 + 5e8 - 8e8) / 3.0,
                                             (-3e8 + 6e8 + 9e8) / 3.0 }));

    EXPECT_FALSE (tr.point_lie_in_plane_tr ({ 1.0, 1.0, 1.0 }));
}

TEST (test_triangle, hard_small)
{
    point_t t1 (1e-9, 2e-9, -3e-9);
    point_t t2 (-4e-9, 5e-9, 6e-9);
    point_t t3 (7e-9, -8e-9, 9e-9);
    triangle_t tr (t1, t2, t3);

    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 1e-9, 2e-9, -3e-9 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ -4e-9, 5e-9, 6e-9 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 7e-9, -8e-9, 9e-9 }));

    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ (1e-9 - 4e-9 + 7e-9) / 3.0,
                                             (2e-9 + 5e-9 - 8e-9) / 3.0,
                                             (-3e-9 + 6e-9 + 9e-9) / 3.0 }));

    EXPECT_FALSE (tr.point_lie_in_plane_tr ({ 1e-9, 2e-9, 1.0 }));
}

TEST (test_triangle, hard_3)
{
    point_t t1 (-1.1, -2.2, -3.3);
    point_t t2 (-4.4, -5.5, -6.6);
    point_t t3 (-7.7, -8.8, -9.9);
    triangle_t tr (t1, t2, t3);

    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ -1.1, -2.2, -3.3 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ -4.4, -5.5, -6.6 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ -7.7, -8.8, -9.9 }));

    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ (-1.1 - 4.4 - 7.7) / 3.0,
                                             (-2.2 - 5.5 - 8.8) / 3.0,
                                             (-3.3 - 6.6 - 9.9) / 3.0 }));

    EXPECT_FALSE (tr.point_lie_in_plane_tr ({ 0.0, 0.0, 0.0 }));
}

// ----------------------------------------------------------------------------------

// ------------------------------TESTING_GET_POINT-----------------------------------

TEST (test_triangle, get_point)
{
    point_t t1 (0.0, 0.5, 1.0);
    point_t t2 (-1.0, 1.5, 2.0);
    point_t t3 (2.0, 2.5, -3.0);
    triangle_t tr (t1, t2, t3);
    EXPECT_DOUBLE_EQ (t1.x_, tr.get_a ().x_);
    EXPECT_DOUBLE_EQ (t1.y_, tr.get_a ().y_);
    EXPECT_DOUBLE_EQ (t1.z_, tr.get_a ().z_);

    EXPECT_DOUBLE_EQ (t2.x_, tr.get_b ().x_);
    EXPECT_DOUBLE_EQ (t2.y_, tr.get_b ().y_);
    EXPECT_DOUBLE_EQ (t2.z_, tr.get_b ().z_);

    EXPECT_DOUBLE_EQ (t3.x_, tr.get_c ().x_);
    EXPECT_DOUBLE_EQ (t3.y_, tr.get_c ().y_);
    EXPECT_DOUBLE_EQ (t3.z_, tr.get_c ().z_);
}

// ----------------------------------------------------------------------------------

// ------------------------------TESTING_INTERSECTION_TR-----------------------------

TEST (intersection_tr, test1)
{
    point_t t1 (0, 0, 1);
    point_t t2 (2, 0, 0);
    point_t t3 (0, 2, -1);

    point_t t4 (0.1, 0.5, -1);
    point_t t5 (0.5, 1, 1);
    point_t t6 (2, 2, 0);

    EXPECT_TRUE (check_intersection_tr_of_line (triangle_t { t1, t2, t3 },
                                                triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test2)
{
    point_t t1 (0, 0, 0);
    point_t t2 (2, 0, 0);
    point_t t3 (0, 2, 0);

    point_t t4 (1, 1, -1);
    point_t t5 (1, 1, 1);
    point_t t6 (2, 2, 0);

    EXPECT_TRUE (check_intersection_tr_of_line (triangle_t { t1, t2, t3 },
                                                triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test3)
{
    point_t t1 (0, 0, 0);
    point_t t2 (1, 0, 0);
    point_t t3 (0, 1, 0);

    point_t t4 (1, 0, 0);
    point_t t5 (0, 1, 0);
    point_t t6 (1, 1, 1);

    EXPECT_TRUE (check_intersection_tr_of_line (triangle_t { t1, t2, t3 },
                                                triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test4)
{
    point_t t1 (0, 0, 0);
    point_t t2 (1, 0, 0);
    point_t t3 (0, 1, 0);

    point_t t4 (0.9, -0.5, 0);
    point_t t5 (-1, 1, 0);
    point_t t6 (1, 1, 1);

    EXPECT_TRUE (check_intersection_tr_of_line (triangle_t { t1, t2, t3 },
                                                triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test5)
{
    point_t t1 (0, 0, 0);
    point_t t2 (1, 0, 0);
    point_t t3 (0, 1, 0);

    point_t t4 (0.4, 0.4, 0);
    point_t t5 (-1, 1, 2);
    point_t t6 (1, 1, 1);

    EXPECT_TRUE (check_intersection_tr_of_line (triangle_t { t1, t2, t3 },
                                                triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test6)
{
    point_t t1 (0, 0, 0);
    point_t t2 (1, 0, 0);
    point_t t3 (0, 1, 0);

    point_t t4 (0.4, 0.4, 0);
    point_t t5 (-1, 1, 0);
    point_t t6 (1, 1, 0);

    EXPECT_TRUE (check_intersection (triangle_t { t1, t2, t3 },
                                     triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test7)
{
    point_t t1 (0, 0, 0);
    point_t t2 (2, 0, 0);
    point_t t3 (0, 2, 0);

    point_t t4 (1, 1, -1);
    point_t t5 (1, 1, 1);
    point_t t6 (2, 2, 0);

    EXPECT_TRUE (check_intersection (triangle_t { t1, t2, t3 },
                                     triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test8)
{
    point_t t1 (0, 0, 0);
    point_t t2 (3, 0, 0);
    point_t t3 (0, 3, 0);

    point_t t4 (1, 1, -1);
    point_t t5 (1, -1, 1);
    point_t t6 (2, 2, 1);

    EXPECT_TRUE (check_intersection (triangle_t { t1, t2, t3 },
                                     triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test9)
{
    point_t t1 (0, 0, 0);
    point_t t2 (1, 0, 0);
    point_t t3 (0, 1, 0);

    point_t t4 (2, 2, 1);
    point_t t5 (3, 2, 1);
    point_t t6 (2, 3, 1);

    EXPECT_FALSE (check_intersection (triangle_t { t1, t2, t3 },
                                      triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test10)
{
    point_t t1 (-2, 0, 0);
    point_t t2 (-1, 1, 0);
    point_t t3 (-1, -1, 0);

    point_t t4 (2, 0, 1);
    point_t t5 (3, 1, 1);
    point_t t6 (3, -1, 1);

    EXPECT_FALSE (check_intersection (triangle_t { t1, t2, t3 },
                                      triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test11)
{
    point_t t1 (0, 0, 0);
    point_t t2 (1, 0, 0);
    point_t t3 (0, 1, 0);

    point_t t4 (0, 0, 0);
    point_t t5 (0, 0, 1);
    point_t t6 (0, 1, 0);

    EXPECT_TRUE (check_intersection (triangle_t { t1, t2, t3 },
                                     triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test12)
{
    point_t t1 (0, 0, 0);
    point_t t2 (1, 0, 0);
    point_t t3 (0, 1, 0);

    point_t t4 (0, 0, 0);
    point_t t5 (1, 0, 0);
    point_t t6 (0, 0, 1);

    EXPECT_TRUE (check_intersection (triangle_t { t1, t2, t3 },
                                     triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test13)
{
    point_t t1 (-1, -1, 0);
    point_t t2 (1, -1, 0);
    point_t t3 (0, 1, 0);

    point_t t4 (0, 0, -1);
    point_t t5 (0, 0, 1);
    point_t t6 (1, 1, 0);

    EXPECT_TRUE (check_intersection (triangle_t { t1, t2, t3 },
                                     triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test14)
{
    point_t t1 (0, 0, 0);
    point_t t2 (1, 0, 0);
    point_t t3 (0, 1, 0);

    point_t t4 (0, 0, 1);
    point_t t5 (1, 0, 1);
    point_t t6 (0, 1, 1);

    EXPECT_FALSE (check_intersection (triangle_t { t1, t2, t3 },
                                      triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test15)
{
    point_t t1 (-2, -2, 0);
    point_t t2 (-1, -2, 0);
    point_t t3 (-2, -1, 0);

    point_t t4 (2, 2, 1);
    point_t t5 (3, 2, -1);
    point_t t6 (2, 3, 1);

    EXPECT_FALSE (check_intersection (triangle_t { t1, t2, t3 },
                                      triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test16)
{
    point_t t1 (0, 0, 0);
    point_t t2 (1, 0, 0);
    point_t t3 (0, 1, 0);

    point_t t4 (2, 2, 1);
    point_t t5 (3, 2, 0);
    point_t t6 (2, 3, 0);

    EXPECT_FALSE (check_intersection (triangle_t { t1, t2, t3 },
                                      triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test17)
{
    point_t t1 (0, 0, 0);
    point_t t2 (1, 0, 0);
    point_t t3 (0, 1, 0);

    point_t t4 (0, 0, 0);
    point_t t5 (-1, 0, 1);
    point_t t6 (0, -1, 1);

    EXPECT_TRUE (check_intersection (triangle_t { t1, t2, t3 },
                                     triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test18)
{
    point_t t1 (0, 0, 0);
    point_t t2 (2, 0, 0.001);
    point_t t3 (0, 2, -0.001);

    point_t t4 (1, 1, -1);
    point_t t5 (1, 1, 1);
    point_t t6 (2, 2, 0);

    EXPECT_TRUE (check_intersection (triangle_t { t1, t2, t3 },
                                     triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test19)
{
    point_t t1 (0, 0, 0);
    point_t t2 (2, 0, 0);
    point_t t3 (0, 2, 0);

    point_t t4 (1, 0, 0);
    point_t t5 (0, 1, 0);
    point_t t6 (2, 2, 0);

    EXPECT_TRUE (check_intersection (triangle_t { t1, t2, t3 },
                                     triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test20)
{
    point_t t1 (0, 0, 0);
    point_t t2 (1, 0, 0);
    point_t t3 (0, 1, 0);

    point_t t4 (2, 2, 0);
    point_t t5 (3, 2, 0);
    point_t t6 (2, 3, 0);

    EXPECT_FALSE (check_intersection (triangle_t { t1, t2, t3 },
                                      triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test21)
{
    point_t t1 (0, 0, 0);
    point_t t2 (1, 0, 0);
    point_t t3 (0, 1, 0);

    point_t t4 (0, 1, 0);
    point_t t5 (1, 0, 0);
    point_t t6 (1, 1, 0);

    EXPECT_TRUE (check_intersection (triangle_t { t1, t2, t3 },
                                     triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test22)
{
    point_t t1 (0, 0, 0);
    point_t t2 (2, 0, 0);
    point_t t3 (0, 2, 0);

    point_t t4 (0.5, 0.5, 0);
    point_t t5 (0.5, 0.5, 1);
    point_t t6 (0.5, 1, -1);

    EXPECT_TRUE (check_intersection (triangle_t { t1, t2, t3 },
                                     triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test23)
{
    point_t t1 (0, 0, 0);
    point_t t2 (2, 0, 0);
    point_t t3 (0, 2, 0);

    point_t t4 (3, 3, -1);
    point_t t5 (4, 3, 1);
    point_t t6 (3, 4, 0);

    EXPECT_FALSE (check_intersection (triangle_t { t1, t2, t3 },
                                      triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test24)
{
    point_t t1 (0, 0, 0);
    point_t t2 (1, 0, 0);
    point_t t3 (0, 1, 0);

    point_t t4 (0.5, 0.5, 0.001);
    point_t t5 (1.5, 0.5, 0.001);
    point_t t6 (0.5, 1.5, 0.001);

    EXPECT_FALSE (check_intersection (triangle_t { t1, t2, t3 },
                                      triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test25)
{
    point_t t1 (0, 0, 0);
    point_t t2 (2, 0, 0);
    point_t t3 (0, 2, 0);

    point_t t4 (1, -1, 0);
    point_t t5 (1, 1, 0);
    point_t t6 (2, 1, 0);

    EXPECT_TRUE (check_intersection (triangle_t { t1, t2, t3 },
                                     triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test26)
{
    point_t t1 (0, 0, 0);
    point_t t2 (1, 0, 0);
    point_t t3 (0, 1, 0);

    point_t t4 (1, 0, 0);
    point_t t5 (1, 1, 1);
    point_t t6 (2, 0, 1);

    EXPECT_TRUE (check_intersection (triangle_t { t1, t2, t3 },
                                     triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test27)
{
    point_t t1 (0, 0, 0);
    point_t t2 (5, 0, 0);
    point_t t3 (0, 5, 0);

    point_t t4 (1, 1, 0);
    point_t t5 (2, 1, 0);
    point_t t6 (1, 2, 0);

    EXPECT_TRUE (check_intersection (triangle_t { t1, t2, t3 },
                                     triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test28)
{
    point_t t1 (0, 0, 0);
    point_t t2 (5, 0, 0);
    point_t t3 (0, 5, 0);

    point_t t4 (1, 1, 0);
    point_t t5 (1, 1, 0);
    point_t t6 (1, 1, 1);

    EXPECT_TRUE (check_intersection (triangle_t { t1, t2, t3 },
                                     triangle_t { t4, t5, t6 }));
}

TEST(intersection_tr, test29)
{
    point_t r1(0, 0, 0);
    point_t r2(5, 0, 0);
    point_t r3(0, 5, 0);

    point_t b1(1, 1, 0); 
    point_t b2(1, 1, 2); 
    point_t b3(2, 2, -2);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test30)
{
    point_t r1(0, 0, 0);
    point_t r2(5, 0, 0);
    point_t r3(0, 5, 0);

    point_t b1(1, 1, 2); 
    point_t b2(2.51, 2.51, 0); 
    point_t b3(2, 2, -2);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test31)
{
    point_t r1(0, 0, 0);
    point_t r2(5, 0, 0);
    point_t r3(0, 5, 0);

    point_t b1(2, 0, 0);
    point_t b2(2, 0, 3);
    point_t b3(2, -10, -3);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test32)
{
    point_t r1(0, 0, 0);
    point_t r2(5, 0, 0);
    point_t r3(0, 5, 0);

    point_t b1(2, -0.1, 0);
    point_t b2(2, 0, 3);
    point_t b3(2, -10, -3);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test33)
{
    point_t r1(0, 0, 0);
    point_t r2(5, 0, 0);
    point_t r3(0, 5, 0);

    point_t b1(0, -0.1, 0);
    point_t b2(2, 0, 3);
    point_t b3(-10, -10, -3);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test34)
{
    point_t r1(0, 0, 0);
    point_t r2(5, 0, 0);
    point_t r3(0, 5, 0);

    point_t b1(0, 0, 0);
    point_t b2(2, 0, 3);
    point_t b3(-10, -10, -3);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test35)
{
    point_t r1(0, 0, 0);
    point_t r2(5, 0, 0);
    point_t r3(0, 5, 0);

    point_t b1(-0.1, 0.1, 0);
    point_t b2(2, 0, 3);
    point_t b3(10, 10, -3);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test36)
{
    point_t r1(0, 0, 0);
    point_t r2(5, 0, 0);
    point_t r3(0, 5, 0);

    point_t b1(2.51, 2.51, 0); 
    point_t b2(1, 1, 2); 
    point_t b3(10, 10, -2);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test37)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b1(0, 2, 3); 
    point_t b2(-0.1, 0.5, 0); 
    point_t b3(2, 0, 2);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test38)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b1(0, 2, 3); 
    point_t b3(-0.1, 0.5, 0); 
    point_t b2(2, 0, 2);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test39)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b2(0, 2, 3); 
    point_t b3(-0.1, 0.5, 0); 
    point_t b1(2, 0, 2);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test40)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b2(0, 2, 3); 
    point_t b1(-0.1, 0.5, 0); 
    point_t b3(2, 0, 2);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test41)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b3(0, 2, 3); 
    point_t b1(-0.1, 0.5, 0); 
    point_t b2(2, 0, 2);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test42)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b3(0, 2, 3); 
    point_t b2(-0.1, 0.5, 0); 
    point_t b1(2, 0, 2);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test43)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b2(0, 2, 3); 
    point_t b3(-0.1, 0.5, 0); 
    point_t b1(2, 0, 2);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test44)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b1(0, 2, 3); 
    point_t b2(0.5, 0.5, 0); 
    point_t b3(2, 0, 2);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test45)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b1(0, 2, 3); 
    point_t b3(0.5, 0.5, 0); 
    point_t b2(2, 0, 2);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test46)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b2(0, 2, 3); 
    point_t b1(0.5, 0.5, 0); 
    point_t b3(2, 0, 2);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test47)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b2(0, 2, 3); 
    point_t b3(0.5, 0.5, 0); 
    point_t b1(2, 0, 2);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test48)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b3(0, 2, 3); 
    point_t b2(0.5, 0.5, 0); 
    point_t b1(2, 0, 2);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test49)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b3(0, 2, 3); 
    point_t b1(0.5, 0.5, 0); 
    point_t b2(2, 0, 2);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test50)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b3(0, 1, -3); 
    point_t b1(3, 3, 0); 
    point_t b2(1, 0, 2);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test51)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b3(0, 1, -3); 
    point_t b2(3, 3, 0); 
    point_t b1(1, 0, 2);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test52)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b1(0, 1, -3); 
    point_t b2(3, 3, 0); 
    point_t b3(1, 0, 2);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test53)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b1(0, 1, -3); 
    point_t b3(3, 3, 0); 
    point_t b2(1, 0, 2);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test54)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b2(0, 1, -3); 
    point_t b1(3, 3, 0); 
    point_t b3(1, 0, 2);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test55)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b2(0, 1, -3); 
    point_t b3(3, 3, 0); 
    point_t b1(1, 0, 2);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test56)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b2(0, 1, -3); 
    point_t b3(3, 3, 0); 
    point_t b1(3, 0, 2);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test57)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b2(0, 1, -3); 
    point_t b1(3, 3, 0); 
    point_t b3(3, 0, 2);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test58)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b1(0, 1, -3); 
    point_t b2(3, 3, 0); 
    point_t b3(3, 0, 2);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test59)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b1(0, 1, -3); 
    point_t b3(3, 3, 0); 
    point_t b2(3, 0, 2);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test60)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b3(0, 1, -3); 
    point_t b2(3, 3, 0); 
    point_t b1(3, 0, 2);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr, test61)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b3(0, 1, -3); 
    point_t b1(3, 3, 0); 
    point_t b2(3, 0, 2);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

// ----------------------------------------------------------------------------------

TEST (point_triangle, test1)
{
    point_t p1 (0, 0, 0);
    point_t p2 (5, 0, 0);
    point_t p3 (0, 5, 0);

    point_t p (1, 1, 0);

    EXPECT_TRUE (check_triangle_point (triangle_t { p1, p2, p3 }, p));
}

TEST (point_triangle, test2)
{
    point_t p1 (0, 0, 0);
    point_t p2 (5, 0, 0);
    point_t p3 (0, 5, 0);

    point_t p (5.001, 0, 0);

    EXPECT_FALSE (check_triangle_point (triangle_t { p1, p2, p3 }, p));
}

TEST (point_triangle, test3)
{
    point_t p1 (-5, 3, 2);
    point_t p2 (11, -6, 10);
    point_t p3 (2, 3, 4);

    point_t p4 (2, 3, 4);
    point_t p5 (11, -6, 10);
    point_t p6 (-5, 3, 2);

    EXPECT_TRUE (check_triangle_point (triangle_t { p1, p2, p3 }, p4));
    EXPECT_TRUE (check_triangle_point (triangle_t { p1, p2, p3 }, p5));
    EXPECT_TRUE (check_triangle_point (triangle_t { p1, p2, p3 }, p6));
}

// ----------------------------------------------------------------------------------

TEST (line_triangle, test1)
{
    point_t p1 (5, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (0, 0, 0);

    point_t p4 (5, 0, 0);
    point_t p5 (10, 10, 10);   

    EXPECT_TRUE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test2)
{
    point_t p1 (5, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (0, 0, 0);

    point_t p4 (1, 1, 0);
    point_t p5 (1, 1, 10);   

    EXPECT_TRUE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test3)
{
    point_t p1 (5, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (0, 0, 0);

    point_t p4 (-1, -1, 0);
    point_t p5 (1, 1, 10);   

    EXPECT_FALSE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test4)
{
    point_t p1 (5, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (0, 0, 0);

    point_t p4 (1, 1, -10);
    point_t p5 (1, 1, 10);   

    EXPECT_TRUE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test5)
{
    point_t p1 (5, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (0, 0, 0);

    point_t p4 (3, 2, 0);
    point_t p5 (1, 1, 10);   

    EXPECT_TRUE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test6)
{
    point_t p1 (5, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (0, 0, 0);

    point_t p4 (3.05, 2, 0);
    point_t p5 (1, 1, 10);   

    EXPECT_FALSE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test7)
{
    point_t p1 (5, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (0, 0, 0);

    point_t p4 (2, 2, -1);
    point_t p5 (1, 1, 10);   

    EXPECT_TRUE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test8)
{
    point_t p1 (5, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (0, 0, 0);

    point_t p4 (0, 0, 10);
    point_t p5 (5, 5, -10);   

    EXPECT_TRUE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test9)
{
    point_t p1 (5, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (0, 0, 0);

    point_t p4 (0, 0, 10.1);
    point_t p5 (5, 5, -10);   

    EXPECT_FALSE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test10)
{
    point_t p1 (4, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (1, 1, 0);

    point_t p4 (3, -1, 0);
    point_t p5 (6, 2, 0);   

    EXPECT_TRUE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test11)
{
    point_t p1 (4, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (1, 1, 0);

    point_t p4 (3, -1.1, 0);
    point_t p5 (6, 2, 0);   

    EXPECT_FALSE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test12)
{
    point_t p1 (4, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (1, 1, 0);

    point_t p4 (1, 1, 0);
    point_t p5 (4, 0, 0);   

    EXPECT_TRUE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test13)
{
    point_t p1 (4, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (1, 1, 0);

    point_t p4 (-2, 2, 0);
    point_t p5 (4, 0, 0);   

    EXPECT_TRUE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test14)
{
    point_t p1 (4, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (1, 1, 0);

    point_t p4 (-2, 2, 0);
    point_t p5 (4, -0.1, 0);   

    EXPECT_FALSE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test15)
{
    point_t p1 (4, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (1, 1, 0);

    point_t p4 (10, 10, 0);
    point_t p5 (0.1, 5, 0);   

    EXPECT_FALSE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test16)
{
    point_t p1 (4, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (1, 1, 0);

    point_t p4 (10, 10, 0);
    point_t p5 (0, 5, 0);   

    EXPECT_TRUE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test17)
{
    point_t p1 (4, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (1, 1, 0);

    point_t p4 (3, 1, 0);
    point_t p5 (2, 2, 0);   

    EXPECT_TRUE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test18)
{
    point_t p1 (4, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (1, 1, 0);

    point_t p4 (0, 1, 0);
    point_t p5 (2, 2, 0);   

    EXPECT_TRUE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test19)
{
    point_t p1 (4, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (1, 1, 0);

    point_t p4 (0, 1, 0);
    point_t p5 (2, 0.6, 0);   

    EXPECT_FALSE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test20)
{
    point_t p1 (4, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (1, 1, 0);

    point_t p4 (0, 0, 0);
    point_t p5 (10, 10, 0);   

    EXPECT_TRUE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test21)
{
    point_t p1 (4, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (1, 1, 0);

    point_t p4 (0.5, 3, 0);
    point_t p5 (0, 0, 0);   

    EXPECT_TRUE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

TEST (line_triangle, test22)
{
    point_t p1 (4, 0, 0);
    point_t p2 (0, 5, 0);
    point_t p3 (1, 1, 0);

    point_t p4 (0.49, 3, 0);
    point_t p5 (0, 0, 0);   

    EXPECT_FALSE (check_triangle_line (triangle_t { p1, p2, p3 }, p4, p5));
}

// ----------------------------------------------------------------------------------

TEST (line_point, test1)
{
    point_t p1 (2, 0, 0);
    point_t p2 (0, 2, 0);

    point_t p3 (1, 1, 0);
    EXPECT_TRUE (check_line_point (p1, p2, p3));
}

TEST(line_point, test2)
{
    point_t A(0, 0, 0);
    point_t B(1, 1, 1);
    point_t P(2, 2, 2);

    EXPECT_FALSE(check_line_point(A, B, P));
}

TEST(line_point, test3)
{
    point_t A(1, 2, 3);
    point_t B(4, 5, 6);
    point_t P(1, 2, 3);

    EXPECT_TRUE(check_line_point(A, B, P));
}

TEST(line_point, test4)
{
    point_t A(1, 2, 3);
    point_t B(4, 5, 6);
    point_t P(4, 5, 6);

    EXPECT_TRUE(check_line_point(A, B, P));
}

TEST(line_point, test5)
{
    point_t A(-1000, -2000, -3000);
    point_t B(-500, -1000, -1500);
    point_t P(-750, -1500, -2250);

    EXPECT_TRUE(check_line_point(A, B, P));
}

TEST(line_point, test6)
{
    point_t A(-1000, -2000, -3000);
    point_t B(-500, -1000, -1500);
    point_t P(-1200, -2400, -3600);

    EXPECT_FALSE(check_line_point(A, B, P));
}

TEST(line_point, test7)
{
    point_t A(0, 0, 0);
    point_t B(1, 1, 1);
    point_t P(1.0 + 1e-10, 1.0 + 1e-10, 1.0 + 1e-10);

    EXPECT_TRUE(check_line_point(A, B, P));
}

TEST(line_point, test8)
{
    point_t A(0, 0, 0);
    point_t B(1, 1, 1);
    point_t P(1000, 1000, 1000);

    EXPECT_FALSE(check_line_point(A, B, P));
}

// ----------------------------------------------------------------------------------

TEST (point_point, test1)
{
    point_t p1 (10, 11, 101);
    point_t p2 (9.99, 11, 101);
    point_t p3 (9.99, 11, 101);

    EXPECT_TRUE (check_point_point (p2, p3));
    EXPECT_FALSE (check_point_point (p1, p2));
}

// ----------------------------------------------------------------------------------

TEST (line_line, test1)
{
    point_t p1 (0, 0, 0);
    point_t p2 (10, 10, 0);
    point_t p3 (10, 0, 0);
    point_t p4 (0, 10, 0);

    EXPECT_TRUE (check_line_line (p1, p2, p3, p4));
}

TEST(line_line, test2) 
{
    point_t A(0, 0, 0);
    point_t B(10, 10, 10);
    point_t C(0, 10, 10);
    point_t D(10, 0, 0);

    EXPECT_TRUE(check_line_line(A, B, C, D));
}

TEST(line_line, test3) 
{
    point_t A(0, 0, 0);
    point_t B(1, 0, 0);
    point_t C(0, 1, 1);
    point_t D(1, 1, 1);

    EXPECT_FALSE(check_line_line(A, B, C, D));
}

TEST(line_line, test4) 
{
    point_t A(0, 0, 0);
    point_t B(5, 5, 5);
    point_t C(3, 3, 3);
    point_t D(8, 8, 8);

    EXPECT_TRUE(check_line_line(A, B, C, D));
}

TEST(line_line, test5) 
{
    point_t A(0, 0, 0);
    point_t B(2, 2, 2);
    point_t C(3, 3, 3);
    point_t D(5, 5, 5);

    EXPECT_FALSE(check_line_line(A, B, C, D));
}

TEST(line_line, test6) 
{
    point_t A(0, 0, 0);
    point_t B(1, 1, 1);
    point_t C(0, 0, 1);
    point_t D(1, 1, 2);

    EXPECT_FALSE(check_line_line(A, B, C, D));
}

TEST (line_line, test7)
{
    point_t p1 (0, 0, 0);
    point_t p2 (10, 10, 10);
    point_t p3 (0, 0, 0);
    point_t p4 (-100, 2, 5);

    EXPECT_TRUE (check_line_line (p1, p2, p3, p4));
}

TEST(line_line, test8) 
{
    point_t A(1.2, 3.5, -2.1);
    point_t B(7.8, 9.6, 4.3);
    point_t C(8.0, 2.1, 0.5);
    point_t D(0.5, 10.0, 3.0);

    EXPECT_FALSE(check_line_line(A, B, C, D));
}

TEST(line_line, test9) 
{
    point_t A(-2.3, 1.1, 0.0);
    point_t B(4.7, 3.9, 2.5);
    point_t C(4.7, 3.9, 2.5);
    point_t D(10.0, -1.2, 3.4);

    EXPECT_TRUE(check_line_line(A, B, C, D));
}

TEST(line_line, test10) 
{
    point_t A(-1.1, 0.0, 2.3);
    point_t B(1.5, 2.7, -0.8);
    point_t C(0.5, 1.3, 1.0);
    point_t D(3.1, 0.9, 0.5);

    EXPECT_FALSE(check_line_line(A, B, C, D));
}

TEST(line_line, test11) 
{
    point_t A(-3.2, 4.1, -1.0);
    point_t B(2.0, 1.5, 3.6);
    point_t C(5.1, 0.3, 2.2);
    point_t D(7.0, -1.0, 4.5);

    EXPECT_FALSE(check_line_line(A, B, C, D));
}

TEST(line_line, test12) 
{
    point_t A(0.7, -2.3, 1.1);
    point_t B(6.5, 3.2, -0.8);
    point_t C(5.0, 0.0, 2.5);
    point_t D(1.2, 4.1, -1.2);

    EXPECT_FALSE(check_line_line(A, B, C, D));
}

// ----------------------------------------------------------------------------------

TEST(intersection_tr1, test1)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b1(1.51, 0.5, 0);
    point_t b2(3, 1, 0);
    point_t b3(1, 5, 0);
    point_t b4(1.5, 0.5, 0);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                    triangle_t{b1, b2, b3}));
    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                    triangle_t{b4, b2, b3}));
}

TEST(intersection_tr1, test2)
{
    point_t r1(2, 0, 0);
    point_t r2(0, 0, 0);
    point_t r3(0, 2, 0);

    point_t b1(2, 0.1, 0);
    point_t b2(3, 1, 0);
    point_t b3(1, 5, 0);
    point_t b4(2, 0, 0);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                    triangle_t{b1, b2, b3}));
    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                    triangle_t{b4, b2, b3}));
}

TEST(intersection_tr1, test3)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b1(2, 0, 0);
    point_t b2(3, 3, 0);
    point_t b3(5, 2, 0);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr1, test4)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b1(0, 0, 0);
    point_t b2(2, 0, 0);
    point_t b3(2, 2, 0);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr1, test5)
{
    point_t r1(0, 0, 0);
    point_t r2(3, 0, 0);
    point_t r3(0, 3, 0);

    point_t b1(2, 0, 0);
    point_t b2(5, 0, 0);
    point_t b3(2, 3, 0);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr1, test6)
{
    point_t r1(0, 0, 0);
    point_t r2(4, 0, 0);
    point_t r3(0, 4, 0);

    point_t b1(1, 1, 0);
    point_t b2(3, 1, 0);
    point_t b3(1, 3, 0);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr1, test7)
{
    point_t r1(0, 0, 0);
    point_t r2(3, 0, 0);
    point_t r3(0, 3, 0);

    point_t b1(1, -1, 0);
    point_t b2(4, 1, 0);
    point_t b3(1, 3, 0);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr1, test8) 
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b1(1.6, 0.5, 0);
    point_t b2(3, 0, 0);
    point_t b3(3, 2, 0);

    EXPECT_FALSE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr1, test9)
{
    point_t r1(0, 0, 0);
    point_t r2(2, 0, 0);
    point_t r3(0, 2, 0);

    point_t b1(-1, -1, 0);
    point_t b2(4, -1, 0);
    point_t b3(0, 5, 0);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

TEST(intersection_tr1, test10)
{
    point_t r1(0, 0, 0);
    point_t r2(3, 0, 0);
    point_t r3(0, 3, 0);

    point_t b1(2, -1, 0);
    point_t b2(4, 2, 0);
    point_t b3(1, 3, 0);

    EXPECT_TRUE(check_intersection(triangle_t{r1, r2, r3},
                                   triangle_t{b1, b2, b3}));
}

// ----------------------------------------------------------------------------------
