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
    EXPECT_DOUBLE_EQ (answer.cor_x, -3.4);
    EXPECT_DOUBLE_EQ (answer.cor_y, 1.6);
    EXPECT_DOUBLE_EQ (answer.cor_z, 0.5);
}

TEST (test_vector, cross_product_zero)
{
    vector_t vec_1 { 0.0, 0.0, 0.0 };
    vector_t vec_2 { 3.0, -5.0, 10.0 };
    vector_t answer = vec_1.cross_product (vec_2);
    EXPECT_DOUBLE_EQ (answer.cor_x, 0.0);
    EXPECT_DOUBLE_EQ (answer.cor_y, 0.0);
    EXPECT_DOUBLE_EQ (answer.cor_z, 0.0);
}

TEST (test_vector, cross_product_parallel)
{
    vector_t vec_1 { 2.0, 2.0, 2.0 };
    vector_t vec_2 { 4.0, 4.0, 4.0 };
    vector_t answer = vec_1.cross_product (vec_2);
    EXPECT_DOUBLE_EQ (answer.cor_x, 0.0);
    EXPECT_DOUBLE_EQ (answer.cor_y, 0.0);
    EXPECT_DOUBLE_EQ (answer.cor_z, 0.0);
}

TEST (test_vector, cross_product_orthogonal)
{
    vector_t vec_1 { 1.0, 0.0, 0.0 };
    vector_t vec_2 { 0.0, 1.0, 0.0 };
    vector_t answer = vec_1.cross_product (vec_2);
    EXPECT_DOUBLE_EQ (answer.cor_x, 0.0);
    EXPECT_DOUBLE_EQ (answer.cor_y, 0.0);
    EXPECT_DOUBLE_EQ (answer.cor_z, 1.0);
}

TEST (test_vector, cross_product_middle)
{
    vector_t vec_1 { 1.25, -3.5, 4.75 };
    vector_t vec_2 { -2.0, 0.5, 3.2 };
    vector_t answer = vec_1.cross_product (vec_2);
    EXPECT_DOUBLE_EQ (answer.cor_x, -13.575);
    EXPECT_DOUBLE_EQ (answer.cor_y, -13.5);
    EXPECT_DOUBLE_EQ (answer.cor_z, -6.375);
}

TEST (test_vector, cross_product_large)
{
    vector_t vec_1 { 1e6, 2e6, -3e6 };
    vector_t vec_2 { -4e6, 5e6, 6e6 };
    vector_t answer = vec_1.cross_product (vec_2);
    EXPECT_DOUBLE_EQ (answer.cor_x, 2.7e13);
    EXPECT_DOUBLE_EQ (answer.cor_y, 6e12);
    EXPECT_DOUBLE_EQ (answer.cor_z, 1.3e13);
}

TEST (test_vector, cross_product_small)
{
    vector_t vec_1 { 0.001, 0.002, 0.003 };
    vector_t vec_2 { 0.004, 0.005, 0.006 };
    vector_t answer = vec_1.cross_product (vec_2);
    EXPECT_DOUBLE_EQ (answer.cor_x, -0.000003);
    EXPECT_DOUBLE_EQ (answer.cor_y, 0.000006);
    EXPECT_DOUBLE_EQ (answer.cor_z, -0.000003);
}

TEST (test_vector, cross_product_negative)
{
    vector_t vec_1 { -1.1, -2.2, -3.3 };
    vector_t vec_2 { -4.4, -5.5, -6.6 };
    vector_t answer = vec_1.cross_product (vec_2);
    EXPECT_DOUBLE_EQ (answer.cor_x, -3.63);
    EXPECT_DOUBLE_EQ (answer.cor_y, 7.26);
    EXPECT_DOUBLE_EQ (answer.cor_z, -3.63);
}

// ----------------------------------------------------------------------------------

// ------------------------------TESTING_POINT_LIE_IN_PLANE_TR-----------------------

TEST (test_triangle, easy_1)
{
    vector_t t1 (5.0, 0.0, 2.7);
    vector_t t2 (1.0, 0.0, 10.3);
    vector_t t3 (2.0, 0.0, 1.0);
    triangle_t tr (t1, t2, t3);
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 13.0, 0.0, 25.9 }));
    EXPECT_FALSE (tr.point_lie_in_plane_tr ({ 1.0, 1.0, 1.0 }));
}

TEST (test_triangle, easy_2)
{
    vector_t t1 (8.9, -7.4, -10.25);
    vector_t t2 (10.0, 0.0, -3.8);
    vector_t t3 (-2.0, 16.0, 17.16);
    triangle_t tr (t1, t2, t3);
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 8.9, -7.4, -10.25 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 10.0, 0.0, -3.8 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ -2.0, 16.0, 17.16 }));
    EXPECT_FALSE (tr.point_lie_in_plane_tr ({ 1.0, 1.0, 1.0 }));
}

TEST (test_triangle, easy_3)
{
    vector_t t1 (0.0, 0.0, 0.0);
    vector_t t2 (1.0, 0.0, 0.0);
    vector_t t3 (0.0, 1.0, 0.0);
    triangle_t tr (t1, t2, t3);

    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 0.0, 0.0, 0.0 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 1.0, 0.0, 0.0 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 0.0, 1.0, 0.0 }));

    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 0.5, 0.5, 0.0 }));

    EXPECT_FALSE (tr.point_lie_in_plane_tr ({ 0.5, 0.5, 1.0 }));
}

TEST (test_triangle, hard_1)
{
    vector_t t1 (1.0, 2.0, 3.0);
    vector_t t2 (-4.0, 5.5, 6.6);
    vector_t t3 (7.7, -8.8, 9.9);
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
    vector_t t1 (-2.5, 4.1, 0.0);
    vector_t t2 (3.3, -1.7, 0.0);
    vector_t t3 (0.0, 0.0, 0.0);
    triangle_t tr (t1, t2, t3);

    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ -2.5, 4.1, 0.0 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 3.3, -1.7, 0.0 }));
    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 0.0, 0.0, 0.0 }));

    EXPECT_TRUE (tr.point_lie_in_plane_tr ({ 100.0, -200.0, 0.0 }));

    EXPECT_FALSE (tr.point_lie_in_plane_tr ({ 0.0, 0.0, 1.0 }));
}

TEST (test_triangle, hard_large)
{
    vector_t t1 (1e8, 2e8, -3e8);
    vector_t t2 (-4e8, 5e8, 6e8);
    vector_t t3 (7e8, -8e8, 9e8);
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
    vector_t t1 (1e-9, 2e-9, -3e-9);
    vector_t t2 (-4e-9, 5e-9, 6e-9);
    vector_t t3 (7e-9, -8e-9, 9e-9);
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
    vector_t t1 (-1.1, -2.2, -3.3);
    vector_t t2 (-4.4, -5.5, -6.6);
    vector_t t3 (-7.7, -8.8, -9.9);
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
    vector_t t1 (0.0, 0.5, 1.0);
    vector_t t2 (-1.0, 1.5, 2.0);
    vector_t t3 (2.0, 2.5, -3.0);
    triangle_t tr (t1, t2, t3);
    EXPECT_DOUBLE_EQ (t1.cor_x, tr.get_a ().cor_x);
    EXPECT_DOUBLE_EQ (t1.cor_y, tr.get_a ().cor_y);
    EXPECT_DOUBLE_EQ (t1.cor_z, tr.get_a ().cor_z);

    EXPECT_DOUBLE_EQ (t2.cor_x, tr.get_b ().cor_x);
    EXPECT_DOUBLE_EQ (t2.cor_y, tr.get_b ().cor_y);
    EXPECT_DOUBLE_EQ (t2.cor_z, tr.get_b ().cor_z);

    EXPECT_DOUBLE_EQ (t3.cor_x, tr.get_c ().cor_x);
    EXPECT_DOUBLE_EQ (t3.cor_y, tr.get_c ().cor_y);
    EXPECT_DOUBLE_EQ (t3.cor_z, tr.get_c ().cor_z);
}

// ----------------------------------------------------------------------------------

// ------------------------------TESTING_INTERSECTION_TR-----------------------------

TEST (intersection_tr, test1)
{
    vector_t t1 (0, 0, 1);
    vector_t t2 (2, 0, 0);
    vector_t t3 (0, 2, -1);

    vector_t t4 (0.1, 0.5, -1);
    vector_t t5 (0.5, 1, 1);
    vector_t t6 (2, 2, 0);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection_tr_of_line(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test2)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (2, 0, 0);
    vector_t t3 (0, 2, 0);

    vector_t t4 (1, 1, -1);
    vector_t t5 (1, 1, 1);
    vector_t t6 (2, 2, 0);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection_tr_of_line(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test3)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (1, 0, 0);
    vector_t t3 (0, 1, 0);

    vector_t t4 (1, 0, 0);
    vector_t t5 (0, 1, 0);
    vector_t t6 (1, 1, 1);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection_tr_of_line(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test4)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (1, 0, 0);
    vector_t t3 (0, 1, 0);

    vector_t t4 (0.9, -0.5, 0);
    vector_t t5 (-1, 1, 0);
    vector_t t6 (1, 1, 1);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection_tr_of_line(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test5)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (1, 0, 0);
    vector_t t3 (0, 1, 0);

    vector_t t4 (0.4, 0.4, 0);
    vector_t t5 (-1, 1, 2);
    vector_t t6 (1, 1, 1);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection_tr_of_line(triangle_t { t4, t5, t6 }));

}

TEST (intersection_tr, test6)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (1, 0, 0);
    vector_t t3 (0, 1, 0);

    vector_t t4 (0.4, 0.4, 0);
    vector_t t5 (-1, 1, 0);
    vector_t t6 (1, 1, 0);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test7)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (2, 0, 0);
    vector_t t3 (0, 2, 0);

    vector_t t4 (1, 1, -1);
    vector_t t5 (1, 1, 1);
    vector_t t6 (2, 2, 0);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test8)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (3, 0, 0);
    vector_t t3 (0, 3, 0);

    vector_t t4 (1, 1, -1);
    vector_t t5 (1, -1, 1);
    vector_t t6 (2, 2, 1);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test9)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (1, 0, 0);
    vector_t t3 (0, 1, 0);

    vector_t t4 (2, 2, 1);
    vector_t t5 (3, 2, 1);
    vector_t t6 (2, 3, 1);

    EXPECT_FALSE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test10)
{
    vector_t t1 (-2, 0, 0);
    vector_t t2 (-1, 1, 0);
    vector_t t3 (-1, -1, 0);

    vector_t t4 (2, 0, 1);
    vector_t t5 (3, 1, 1);
    vector_t t6 (3, -1, 1);

    EXPECT_FALSE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test11)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (1, 0, 0);
    vector_t t3 (0, 1, 0);

    vector_t t4 (0, 0, 0);
    vector_t t5 (0, 0, 1);
    vector_t t6 (0, 1, 0);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test12)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (1, 0, 0);
    vector_t t3 (0, 1, 0);

    vector_t t4 (0, 0, 0);
    vector_t t5 (1, 0, 0);
    vector_t t6 (0, 0, 1);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test13)
{
    vector_t t1 (-1, -1, 0);
    vector_t t2 (1, -1, 0);
    vector_t t3 (0, 1, 0);

    vector_t t4 (0, 0, -1);
    vector_t t5 (0, 0, 1);
    vector_t t6 (1, 1, 0);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test14)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (1, 0, 0);
    vector_t t3 (0, 1, 0);

    vector_t t4 (0, 0, 1);
    vector_t t5 (1, 0, 1);
    vector_t t6 (0, 1, 1);

    EXPECT_FALSE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test15)
{
    vector_t t1 (-2, -2, 0);
    vector_t t2 (-1, -2, 0);
    vector_t t3 (-2, -1, 0);

    vector_t t4 (2, 2, 1);
    vector_t t5 (3, 2, -1);
    vector_t t6 (2, 3, 1);

    EXPECT_FALSE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test16)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (1, 0, 0);
    vector_t t3 (0, 1, 0);

    vector_t t4 (2, 2, 1);
    vector_t t5 (3, 2, 0);
    vector_t t6 (2, 3, 0);

    EXPECT_FALSE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test17)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (1, 0, 0);
    vector_t t3 (0, 1, 0);

    vector_t t4 (0, 0, 0);
    vector_t t5 (-1, 0, 1);
    vector_t t6 (0, -1, 1);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test18)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (2, 0, 0.001);
    vector_t t3 (0, 2, -0.001);

    vector_t t4 (1, 1, -1);
    vector_t t5 (1, 1, 1);
    vector_t t6 (2, 2, 0);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test19)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (2, 0, 0);
    vector_t t3 (0, 2, 0);

    vector_t t4 (1, 0, 0);
    vector_t t5 (0, 1, 0);
    vector_t t6 (2, 2, 0);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test20)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (1, 0, 0);
    vector_t t3 (0, 1, 0);

    vector_t t4 (2, 2, 0);
    vector_t t5 (3, 2, 0);
    vector_t t6 (2, 3, 0);

    EXPECT_FALSE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test21)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (1, 0, 0);
    vector_t t3 (0, 1, 0);

    vector_t t4 (0, 1, 0);
    vector_t t5 (1, 0, 0);
    vector_t t6 (1, 1, 0);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));

}

TEST (intersection_tr, test22)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (2, 0, 0);
    vector_t t3 (0, 2, 0);

    vector_t t4 (0.5, 0.5, 0);
    vector_t t5 (0.5, 0.5, 1);
    vector_t t6 (0.5, 1, -1);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test23)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (2, 0, 0);
    vector_t t3 (0, 2, 0);

    vector_t t4 (3, 3, -1);
    vector_t t5 (4, 3, 1);
    vector_t t6 (3, 4, 0);

    EXPECT_FALSE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test24)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (1, 0, 0);
    vector_t t3 (0, 1, 0);

    vector_t t4 (0.5, 0.5, 0.001);
    vector_t t5 (1.5, 0.5, 0.001);
    vector_t t6 (0.5, 1.5, 0.001);

    EXPECT_FALSE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test25)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (2, 0, 0);
    vector_t t3 (0, 2, 0);

    vector_t t4 (1, -1, 0);
    vector_t t5 (1, 1, 0);
    vector_t t6 (2, 1, 0);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test26)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (1, 0, 0);
    vector_t t3 (0, 1, 0);

    vector_t t4 (1, 0, 0);
    vector_t t5 (1, 1, 1);
    vector_t t6 (2, 0, 1);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test27)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (5, 0, 0);
    vector_t t3 (0, 5, 0);

    vector_t t4 (1, 1, 0);
    vector_t t5 (2, 1, 0);
    vector_t t6 (1, 2, 0);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST (intersection_tr, test28)
{
    vector_t t1 (0, 0, 0);
    vector_t t2 (5, 0, 0);
    vector_t t3 (0, 5, 0);

    vector_t t4 (1, 1, 0);
    vector_t t5 (1, 1, 0);
    vector_t t6 (1, 1, 1);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST(intersection_tr, test29)
{
    vector_t t1(0, 0, 0);
    vector_t t2(5, 0, 0);
    vector_t t3(0, 5, 0);

    vector_t t4(1, 1, 0); 
    vector_t t5(1, 1, 2); 
    vector_t t6(2, 2, -2);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST(intersection_tr, test30)
{
    vector_t t1(0, 0, 0);
    vector_t t2(5, 0, 0);
    vector_t t3(0, 5, 0);

    vector_t t4(1, 1, 2); 
    vector_t t5(2.51, 2.51, 0); 
    vector_t t6(2, 2, -2);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST(intersection_tr, test31)
{
    vector_t t1(0, 0, 0);
    vector_t t2(5, 0, 0);
    vector_t t3(0, 5, 0);

    vector_t t4(2, 0, 0);
    vector_t t5(2, 0, 3);
    vector_t t6(2, -10, -3);

    EXPECT_TRUE ((triangle_t { t1, t2, t3 }).check_intersection(triangle_t { t4, t5, t6 }));
}

TEST(intersection_tr, test32)
{
    vector_t r1(0, 0, 0);
    vector_t r2(5, 0, 0);
    vector_t r3(0, 5, 0);

    vector_t b1(2, -0.1, 0);
    vector_t b2(2, 0, 3);
    vector_t b3(2, -10, -3);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test33)
{
    vector_t r1(0, 0, 0);
    vector_t r2(5, 0, 0);
    vector_t r3(0, 5, 0);

    vector_t b1(0, -0.1, 0);
    vector_t b2(2, 0, 3);
    vector_t b3(-10, -10, -3);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test34)
{
    vector_t r1(0, 0, 0);
    vector_t r2(5, 0, 0);
    vector_t r3(0, 5, 0);

    vector_t b1(0, 0, 0);
    vector_t b2(2, 0, 3);
    vector_t b3(-10, -10, -3);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test35)
{
    vector_t r1(0, 0, 0);
    vector_t r2(5, 0, 0);
    vector_t r3(0, 5, 0);

    vector_t b1(-0.1, 0.1, 0);
    vector_t b2(2, 0, 3);
    vector_t b3(10, 10, -3);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test36)
{
    vector_t r1(0, 0, 0);
    vector_t r2(5, 0, 0);
    vector_t r3(0, 5, 0);

    vector_t b1(2.51, 2.51, 0); 
    vector_t b2(1, 1, 2); 
    vector_t b3(10, 10, -2);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test37)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b1(0, 2, 3); 
    vector_t b2(-0.1, 0.5, 0); 
    vector_t b3(2, 0, 2);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test38)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b1(0, 2, 3); 
    vector_t b3(-0.1, 0.5, 0); 
    vector_t b2(2, 0, 2);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test39)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b2(0, 2, 3); 
    vector_t b3(-0.1, 0.5, 0); 
    vector_t b1(2, 0, 2);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test40)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b2(0, 2, 3); 
    vector_t b1(-0.1, 0.5, 0); 
    vector_t b3(2, 0, 2);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test41)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b3(0, 2, 3); 
    vector_t b1(-0.1, 0.5, 0); 
    vector_t b2(2, 0, 2);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test42)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b3(0, 2, 3); 
    vector_t b2(-0.1, 0.5, 0); 
    vector_t b1(2, 0, 2);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test43)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b2(0, 2, 3); 
    vector_t b3(-0.1, 0.5, 0); 
    vector_t b1(2, 0, 2);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test44)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b1(0, 2, 3); 
    vector_t b2(0.5, 0.5, 0); 
    vector_t b3(2, 0, 2);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test45)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b1(0, 2, 3); 
    vector_t b3(0.5, 0.5, 0); 
    vector_t b2(2, 0, 2);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test46)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b2(0, 2, 3); 
    vector_t b1(0.5, 0.5, 0); 
    vector_t b3(2, 0, 2);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test47)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b2(0, 2, 3); 
    vector_t b3(0.5, 0.5, 0); 
    vector_t b1(2, 0, 2);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test48)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b3(0, 2, 3); 
    vector_t b2(0.5, 0.5, 0); 
    vector_t b1(2, 0, 2);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test49)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b3(0, 2, 3); 
    vector_t b1(0.5, 0.5, 0); 
    vector_t b2(2, 0, 2);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test50)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b3(0, 1, -3); 
    vector_t b1(3, 3, 0); 
    vector_t b2(1, 0, 2);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test51)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b3(0, 1, -3); 
    vector_t b2(3, 3, 0); 
    vector_t b1(1, 0, 2);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test52)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b1(0, 1, -3); 
    vector_t b2(3, 3, 0); 
    vector_t b3(1, 0, 2);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test53)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b1(0, 1, -3); 
    vector_t b3(3, 3, 0); 
    vector_t b2(1, 0, 2);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test54)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b2(0, 1, -3); 
    vector_t b1(3, 3, 0); 
    vector_t b3(1, 0, 2);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test55)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b2(0, 1, -3); 
    vector_t b3(3, 3, 0); 
    vector_t b1(1, 0, 2);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test56)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b2(0, 1, -3); 
    vector_t b3(3, 3, 0); 
    vector_t b1(3, 0, 2);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test57)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b2(0, 1, -3); 
    vector_t b1(3, 3, 0); 
    vector_t b3(3, 0, 2);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test58)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b1(0, 1, -3); 
    vector_t b2(3, 3, 0); 
    vector_t b3(3, 0, 2);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test59)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b1(0, 1, -3); 
    vector_t b3(3, 3, 0); 
    vector_t b2(3, 0, 2);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test60)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b3(0, 1, -3); 
    vector_t b2(3, 3, 0); 
    vector_t b1(3, 0, 2);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr, test61)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b3(0, 1, -3); 
    vector_t b1(3, 3, 0); 
    vector_t b2(3, 0, 2);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

// ----------------------------------------------------------------------------------

TEST (vector_triangle, test1)
{
    vector_t p1 (0, 0, 0);
    vector_t p2 (5, 0, 0);
    vector_t p3 (0, 5, 0);

    vector_t p (1, 1, 0);

    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_point (p));
}

TEST (vector_triangle, test2)
{
    vector_t p1 (0, 0, 0);
    vector_t p2 (5, 0, 0);
    vector_t p3 (0, 5, 0);

    vector_t p (5.001, 0, 0);

    EXPECT_FALSE ((triangle_t { p1, p2, p3 }).check_triangle_point (p));
}

TEST (vector_triangle, test3)
{
    vector_t p1 (-5, 3, 2);
    vector_t p2 (11, -6, 10);
    vector_t p3 (2, 3, 4);

    vector_t p4 (2, 3, 4);
    vector_t p5 (11, -6, 10);
    vector_t p6 (-5, 3, 2);

    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_point (p4));
    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_point (p5));
    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_point (p6));

}

// ----------------------------------------------------------------------------------

TEST (line_triangle, test1)
{
    vector_t p1 (5, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (0, 0, 0);

    vector_t p4 (5, 0, 0);
    vector_t p5 (10, 10, 10);   

    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test2)
{
    vector_t p1 (5, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (0, 0, 0);

    vector_t p4 (1, 1, 0);
    vector_t p5 (1, 1, 10);   

    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test3)
{
    vector_t p1 (5, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (0, 0, 0);

    vector_t p4 (-1, -1, 0);
    vector_t p5 (1, 1, 10);   

    EXPECT_FALSE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test4)
{
    vector_t p1 (5, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (0, 0, 0);

    vector_t p4 (1, 1, -10);
    vector_t p5 (1, 1, 10);   

    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test5)
{
    vector_t p1 (5, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (0, 0, 0);

    vector_t p4 (3, 2, 0);
    vector_t p5 (1, 1, 10);   

    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test6)
{
    vector_t p1 (5, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (0, 0, 0);

    vector_t p4 (3.05, 2, 0);
    vector_t p5 (1, 1, 10);   

    EXPECT_FALSE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test7)
{
    vector_t p1 (5, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (0, 0, 0);

    vector_t p4 (2, 2, -1);
    vector_t p5 (1, 1, 10);   

    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test8)
{
    vector_t p1 (5, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (0, 0, 0);

    vector_t p4 (0, 0, 10);
    vector_t p5 (5, 5, -10);   

    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test9)
{
    vector_t p1 (5, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (0, 0, 0);

    vector_t p4 (0, 0, 10.1);
    vector_t p5 (5, 5, -10);   

    EXPECT_FALSE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test10)
{
    vector_t p1 (4, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (1, 1, 0);

    vector_t p4 (3, -1, 0);
    vector_t p5 (6, 2, 0);   

    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test11)
{
    vector_t p1 (4, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (1, 1, 0);

    vector_t p4 (3, -1.1, 0);
    vector_t p5 (6, 2, 0);   

    EXPECT_FALSE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test12)
{
    vector_t p1 (4, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (1, 1, 0);

    vector_t p4 (1, 1, 0);
    vector_t p5 (4, 0, 0);   

    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test13)
{
    vector_t p1 (4, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (1, 1, 0);

    vector_t p4 (-2, 2, 0);
    vector_t p5 (4, 0, 0);   

    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test14)
{
    vector_t p1 (4, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (1, 1, 0);

    vector_t p4 (-2, 2, 0);
    vector_t p5 (4, -0.1, 0);   

    EXPECT_FALSE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test15)
{
    vector_t p1 (4, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (1, 1, 0);

    vector_t p4 (10, 10, 0);
    vector_t p5 (0.1, 5, 0);   

    EXPECT_FALSE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test16)
{
    vector_t p1 (4, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (1, 1, 0);

    vector_t p4 (10, 10, 0);
    vector_t p5 (0, 5, 0);   

    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test17)
{
    vector_t p1 (4, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (1, 1, 0);

    vector_t p4 (3, 1, 0);
    vector_t p5 (2, 2, 0);   

    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test18)
{
    vector_t p1 (4, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (1, 1, 0);

    vector_t p4 (0, 1, 0);
    vector_t p5 (2, 2, 0);   

    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test19)
{
    vector_t p1 (4, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (1, 1, 0);

    vector_t p4 (0, 1, 0);
    vector_t p5 (2, 0.6, 0);   

    EXPECT_FALSE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test20)
{
    vector_t p1 (4, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (1, 1, 0);

    vector_t p4 (0, 0, 0);
    vector_t p5 (10, 10, 0);   

    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test21)
{
    vector_t p1 (4, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (1, 1, 0);

    vector_t p4 (0.5, 3, 0);
    vector_t p5 (0, 0, 0);   

    EXPECT_TRUE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

TEST (line_triangle, test22)
{
    vector_t p1 (4, 0, 0);
    vector_t p2 (0, 5, 0);
    vector_t p3 (1, 1, 0);

    vector_t p4 (0.49, 3, 0);
    vector_t p5 (0, 0, 0);   

    EXPECT_FALSE ((triangle_t { p1, p2, p3 }).check_triangle_line (p4, p5));
}

// ----------------------------------------------------------------------------------

TEST (line_point, test1)
{
    vector_t p1 (2, 0, 0);
    vector_t p2 (0, 2, 0);

    vector_t p3 (1, 1, 0);
    EXPECT_TRUE (triangle_t::check_line_point (p1, p2, p3));
}

TEST(line_point, test2)
{
    vector_t A(0, 0, 0);
    vector_t B(1, 1, 1);
    vector_t P(2, 2, 2);

    EXPECT_FALSE(triangle_t::check_line_point(A, B, P));
}

TEST(line_point, test3)
{
    vector_t A(1, 2, 3);
    vector_t B(4, 5, 6);
    vector_t P(1, 2, 3);

    EXPECT_TRUE(triangle_t::check_line_point(A, B, P));
}

TEST(line_point, test4)
{
    vector_t A(1, 2, 3);
    vector_t B(4, 5, 6);
    vector_t P(4, 5, 6);

    EXPECT_TRUE(triangle_t::check_line_point(A, B, P));
}

TEST(line_point, test5)
{
    vector_t A(-1000, -2000, -3000);
    vector_t B(-500, -1000, -1500);
    vector_t P(-750, -1500, -2250);

    EXPECT_TRUE(triangle_t::check_line_point(A, B, P));
}

TEST(line_point, test6)
{
    vector_t A(-1000, -2000, -3000);
    vector_t B(-500, -1000, -1500);
    vector_t P(-1200, -2400, -3600);

    EXPECT_FALSE(triangle_t::check_line_point(A, B, P));
}

TEST(line_point, test7)
{
    vector_t A(0, 0, 0);
    vector_t B(1, 1, 1);
    vector_t P(1.0 + 1e-10, 1.0 + 1e-10, 1.0 + 1e-10);

    EXPECT_TRUE(triangle_t::check_line_point(A, B, P));
}

TEST(line_point, test8)
{
    vector_t A(0, 0, 0);
    vector_t B(1, 1, 1);
    vector_t P(1000, 1000, 1000);

    EXPECT_FALSE(triangle_t::check_line_point(A, B, P));
}

// ----------------------------------------------------------------------------------

TEST (point_point, test1)
{
    vector_t p1 (10, 11, 101);
    vector_t p2 (9.99, 11, 101);
    vector_t p3 (9.99, 11, 101);

    EXPECT_TRUE (triangle_t::check_point_point (p2, p3));
    EXPECT_FALSE (triangle_t::check_point_point (p1, p2));
}

// ----------------------------------------------------------------------------------

TEST (line_line, test1)
{
    vector_t p1 (0, 0, 0);
    vector_t p2 (10, 10, 0);
    vector_t p3 (10, 0, 0);
    vector_t p4 (0, 10, 0);

    EXPECT_TRUE (triangle_t::check_line_line (p1, p2, p3, p4));
}

TEST(line_line, test2) 
{
    vector_t A(0, 0, 0);
    vector_t B(10, 10, 10);
    vector_t C(0, 10, 10);
    vector_t D(10, 0, 0);

    EXPECT_TRUE(triangle_t::check_line_line(A, B, C, D));
}

TEST(line_line, test3) 
{
    vector_t A(0, 0, 0);
    vector_t B(1, 0, 0);
    vector_t C(0, 1, 1);
    vector_t D(1, 1, 1);

    EXPECT_FALSE(triangle_t::check_line_line(A, B, C, D));
}

TEST(line_line, test4) 
{
    vector_t A(0, 0, 0);
    vector_t B(5, 5, 5);
    vector_t C(3, 3, 3);
    vector_t D(8, 8, 8);

    EXPECT_TRUE(triangle_t::check_line_line(A, B, C, D));
}

TEST(line_line, test5) 
{
    vector_t A(0, 0, 0);
    vector_t B(2, 2, 2);
    vector_t C(3, 3, 3);
    vector_t D(5, 5, 5);

    EXPECT_FALSE(triangle_t::check_line_line(A, B, C, D));
}

TEST(line_line, test6) 
{
    vector_t A(0, 0, 0);
    vector_t B(1, 1, 1);
    vector_t C(0, 0, 1);
    vector_t D(1, 1, 2);

    EXPECT_FALSE(triangle_t::check_line_line(A, B, C, D));
}

TEST (line_line, test7)
{
    vector_t p1 (0, 0, 0);
    vector_t p2 (10, 10, 10);
    vector_t p3 (0, 0, 0);
    vector_t p4 (-100, 2, 5);

    EXPECT_TRUE (triangle_t::check_line_line (p1, p2, p3, p4));
}

TEST(line_line, test8) 
{
    vector_t A(1.2, 3.5, -2.1);
    vector_t B(7.8, 9.6, 4.3);
    vector_t C(8.0, 2.1, 0.5);
    vector_t D(0.5, 10.0, 3.0);

    EXPECT_FALSE(triangle_t::check_line_line(A, B, C, D));
}

TEST(line_line, test9) 
{
    vector_t A(-2.3, 1.1, 0.0);
    vector_t B(4.7, 3.9, 2.5);
    vector_t C(4.7, 3.9, 2.5);
    vector_t D(10.0, -1.2, 3.4);

    EXPECT_TRUE(triangle_t::check_line_line(A, B, C, D));
}

TEST(line_line, test10) 
{
    vector_t A(-1.1, 0.0, 2.3);
    vector_t B(1.5, 2.7, -0.8);
    vector_t C(0.5, 1.3, 1.0);
    vector_t D(3.1, 0.9, 0.5);

    EXPECT_FALSE(triangle_t::check_line_line(A, B, C, D));
}

TEST(line_line, test11) 
{
    vector_t A(-3.2, 4.1, -1.0);
    vector_t B(2.0, 1.5, 3.6);
    vector_t C(5.1, 0.3, 2.2);
    vector_t D(7.0, -1.0, 4.5);

    EXPECT_FALSE(triangle_t::check_line_line(A, B, C, D));
}

TEST(line_line, test12) 
{
    vector_t A(0.7, -2.3, 1.1);
    vector_t B(6.5, 3.2, -0.8);
    vector_t C(5.0, 0.0, 2.5);
    vector_t D(1.2, 4.1, -1.2);

    EXPECT_FALSE(triangle_t::check_line_line(A, B, C, D));
}

// ----------------------------------------------------------------------------------

TEST(intersection_tr1, test1)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b1(1.51, 0.5, 0);
    vector_t b2(3, 1, 0);
    vector_t b3(1, 5, 0);
    vector_t b4(1.5, 0.5, 0);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b4, b2, b3 }));
}

TEST(intersection_tr1, test2)
{
    vector_t r1(2, 0, 0);
    vector_t r2(0, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b1(2, 0.1, 0);
    vector_t b2(3, 1, 0);
    vector_t b3(1, 5, 0);
    vector_t b4(2, 0, 0);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b4, b2, b3 }));
}

TEST(intersection_tr1, test3)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b1(2, 0, 0);
    vector_t b2(3, 3, 0);
    vector_t b3(5, 2, 0);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr1, test4)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b1(0, 0, 0);
    vector_t b2(2, 0, 0);
    vector_t b3(2, 2, 0);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr1, test5)
{
    vector_t r1(0, 0, 0);
    vector_t r2(3, 0, 0);
    vector_t r3(0, 3, 0);

    vector_t b1(2, 0, 0);
    vector_t b2(5, 0, 0);
    vector_t b3(2, 3, 0);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr1, test6)
{
    vector_t r1(0, 0, 0);
    vector_t r2(4, 0, 0);
    vector_t r3(0, 4, 0);

    vector_t b1(1, 1, 0);
    vector_t b2(3, 1, 0);
    vector_t b3(1, 3, 0);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr1, test7)
{
    vector_t r1(0, 0, 0);
    vector_t r2(3, 0, 0);
    vector_t r3(0, 3, 0);

    vector_t b1(1, -1, 0);
    vector_t b2(4, 1, 0);
    vector_t b3(1, 3, 0);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr1, test8) 
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b1(1.6, 0.5, 0);
    vector_t b2(3, 0, 0);
    vector_t b3(3, 2, 0);

    EXPECT_FALSE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr1, test9)
{
    vector_t r1(0, 0, 0);
    vector_t r2(2, 0, 0);
    vector_t r3(0, 2, 0);

    vector_t b1(-1, -1, 0);
    vector_t b2(4, -1, 0);
    vector_t b3(0, 5, 0);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

TEST(intersection_tr1, test10)
{
    vector_t r1(0, 0, 0);
    vector_t r2(3, 0, 0);
    vector_t r3(0, 3, 0);

    vector_t b1(2, -1, 0);
    vector_t b2(4, 2, 0);
    vector_t b3(1, 3, 0);

    EXPECT_TRUE ((triangle_t { r1, r2, r3 }).check_intersection(triangle_t { b1, b2, b3 }));
}

// ----------------------------------------------------------------------------------
