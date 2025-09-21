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

// ----------------------------------------------------------------------------------
