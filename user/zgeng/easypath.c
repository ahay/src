/* Easy path seislet transform */
/*
  Copyright (C) 2018 University of Texas at Austin
   
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
   
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
   
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>

#include "easypath.h"

static int delete_element(int *arr, int n, int x)
/* delete a given element x from an array */
{
    int i=0, j;
    while (i < n) {
        for (i = 0; i < n; i++)
            if (arr[i] == x)
                break;

        if (i < n) {
            n = n - 1;
            for (j = i; j < n; j++)
                arr[j] = arr[j+1];
        }
    }
    return n;
}

static int delete_elements(int *arr1, int n1, int *arr2, int n2)
/* delete several elements in arr2 from arr1 */
{
    int i;
    for (i = 0; i < n2; i++)
        n1 = delete_element(arr1, n1, arr2[i]);

    return n1;
}

static int find_union(int *array1, int n1, int *array2, int n2, int **union_array)
/* find the union of two sorted array */
{
    int i = 0, j = 0, k = 0;
    int *array3 = sf_intalloc(n1+n2);
    while ((i < n1) && (j < n2)) {
        if (array1[i] < array2[j]) 
            array3[k++] = array1[i++];
        else if (array1[i] > array2[j]) 
            array3[k++] = array2[j++];
        else {
            array3[k++] = array1[i++];
            j++;
        }
    }

    if (i == n1) {
        while (j < n2) 
            array3[k++] = array2[j++];
    } else {
        while (i < n1) 
            array3[k++] = array1[i++];
    }

    *union_array = array3;
    return k;
}

static void clip(int *arr1, int n, int low, int high)
/* find values which are out of range and set them as 0 */
{
    int i;
    for (i = 0; i < n; i++) {
        if (arr1[i] < low)
            arr1[i] = 0;
        else if (arr1[i] > high)
            arr1[i] = 0;
    }
}

// return value is the number of neighbours
static int find_neighbours(int *indexes,           /* one set of indexes */ 
                    int num_indexes,        /* the number of indexes */ 
                    int **neighbours_all,   /* neighbours of all indexes; here pass address of *neighbours_all */
                    int n1, int n2          /* size of the data */)
/* find neighbours of all indexes */
{
    int *neighbours = NULL; // array combining neighbours of each index and pass to **neighbours_all
    int num_neighbours = 0; // return value
    int length;             // the number of neighbours of one index
    int *a;                 // temporary arrary to store the neighbours of one index
    int i;

    // get union of neighbours of each index and put it to neighbours array
    for (i = 0; i < num_indexes; i++) {
        if (indexes[i]%n1 == 0) { // top boundary

            length = 5;
            a = sf_intalloc(length);
            // possible indexes
            int b[] = {indexes[i]-n1-1, indexes[i]-n1, indexes[i]-1,
                indexes[i]+n1-1, indexes[i]+n1};
            memcpy(a, b, sizeof(b));

        } else if (indexes[i]%n1 == 1) { // bottom boundary

            length = 5;
            a = sf_intalloc(length);
            // possible indexes
            int b[] = {indexes[i]-n1, indexes[i]-n1+1, indexes[i]+1,
                indexes[i]+n1, indexes[i]+n1+1};
            memcpy(a, b, sizeof(b));

        } else {

            length = 8;
            a = sf_intalloc(length);
            // possible indexes
            int b[] = {indexes[i]-n1-1, indexes[i]-n1, indexes[i]-n1+1,
                indexes[i]-1, indexes[i]+1, indexes[i]+n1-1, 
                indexes[i]+n1, indexes[i]+n1+1};
            memcpy(a, b, sizeof(b));

        }

        // find union between neighbours of previous indexes and neighbours of the current index
        num_neighbours = find_union(neighbours, num_neighbours, a, length, &neighbours);
        free(a);
    }

    // clip values outside of (0, n1*n2) to 0 
    clip(neighbours, num_neighbours, 0, n1*n2);

    // delete index 0
    num_neighbours = delete_element(neighbours, num_neighbours, 0);

    // delete the same indexes 
    num_neighbours = delete_elements(neighbours, num_neighbours, indexes, num_indexes);
    
    *neighbours_all = neighbours;
    
    return num_neighbours;

}

static void find_indexes(int *arr1, int n1, int x, int *indexes)
/* find all indexes belonging to xth index set */
{
    int i, j = 0;
    for (i = 0; i < n1; i++) {
        if (arr1[i] == x)
            indexes[j++] = i+1;
    }
}

static void used_neighbours(int *arr1, int n1, int *a1)
/* find indexes of used neighbours and set them as 0 */
{
    int i;
    for (i = 0; i < n1; i++) {
        if (a1[arr1[i]-1] != 0)
            arr1[i] = 0;
    }
}

static int find_left_indexes(int *arr1, int n1, int **left)
/* find indexes haven't been selected */
{
    int i, j = 0;
    int *arr2 = sf_intalloc(n1);
    for (i = 0; i < n1; i++) {
        if (arr1[i] == 0)
            arr2[j++] = i+1;
    }

    *left = arr2;

    return j;
}

// arr3 is the index of arr2; arr2 is the index of arr1
/* compute the absolute difference between one value and an array */
/*static float* diff(float *arr1, int x, int *arr2, int *arr3, int n3)
{
    float *arr4 = sf_floatalloc(n3);
    for (int i = 0; i < n3; i++) 
        arr4[i] = fabs(arr1[x-1] - arr1[arr2[arr3[i]-1]-1]);

    return arr4;
}*/

static float corr_coeff(float *arr1, float *arr2, int nt)
/* compute the correlation coefficient of two traces */
{
    int i;
    float simi, square_sum1=0., square_sum2=0., multi_sum=0.;

    for (i = 0; i < nt; i++) {
        square_sum1 += arr1[i]*arr1[i];
        square_sum2 += arr2[i]*arr2[i];
        multi_sum   += arr1[i]*arr2[i];
    }
    simi = multi_sum*multi_sum/(square_sum1*square_sum2 + 1e-20);

    return simi;
}

static float* diff(float **arr1, int x, int *arr2, int *arr3, int n3, int nt)
/* compute the ute difference between one value and an array */
{
    int i;
    float *arr4 = sf_floatalloc(n3);
    for (i = 0; i < n3; i++) 
        arr4[i] = corr_coeff(arr1[x-1], arr1[arr2[arr3[i]-1]-1], nt);

    return arr4;
}

static int find_max(float *arr, int n1)
/* find the location of the maximum in an array */
{
    float max = arr[0];
    int i, loc_max = 0;
    for (i = 0; i < n1; i++) {
        if (arr[i] > max) {
            max = arr[i];
            loc_max = i;
        }
    }

    return loc_max;
}

int *pathway(float **data,      /* data array */
            int *index_sets,    /* index sets; it says which part the data value belongs to */
            int nt,             /* length of time samples */
            int n1, int n2      /* size of the index sets */)
/*< find the path for the seislet transform >*/
{
    int i, j, loc_max;
    int iterations = 1;
    int num_neighbours;
    int *neighbours;
    int n = n1*n2*nt;

    // vector indicating the path
    int *path           = sf_intalloc(n1*n2);     
    // new indexes set; connect two neighbour sets together
    int *index_sets_new = sf_intalloc(n1*n2);     
    // temporary array to store the difference between neighbours and the current index
    float *diff_temp; 
    // number of indexes for each index set
    int num_indexes = (int)pow(2, iterations-1);
    int indexes[num_indexes];

    /* initialize */
    for (i = 0; i < n1*n2; i++)
        path[i] = 0;
    for (i = 0; i < n1*n2; i++)
        index_sets_new[i] = 0;

    /* always start from the first index */
    path[0] = 1;


    for (i = 0; i < n1*n2-1; i++) {
        
        /* find all indexes belonging to path[i]th index set */
        find_indexes(index_sets, n1*n2, path[i], indexes);

        /* connect two neighbour indexes sets together */
        for (j = 0; j < num_indexes; j++)
            index_sets_new[indexes[j]-1] = floor(i/2)+1;

        /* find neighbours around all indexes */
        num_neighbours = find_neighbours(indexes, num_indexes, &neighbours, n1, n2);

        /* don't use neighbours that are already used in the path */
        used_neighbours(neighbours, num_neighbours, index_sets_new);

        /* delete 0 indexes */
        num_neighbours = delete_element(neighbours, num_neighbours, 0);

        /* if there are no neighbour then take all left indexes */
        if (num_neighbours == 0) 
            num_neighbours = find_left_indexes(index_sets_new, n1*n2, &neighbours);

        /* get the difference of function values between neighbours and the current index */
        diff_temp = diff(data, path[i], index_sets, neighbours, num_neighbours, nt);

        /* get the index of the maximum */
        loc_max = find_max(diff_temp, num_neighbours);

        /* find the index corresponding to the maximum and store the data value */
        path[i+1] = index_sets[neighbours[loc_max]-1];
    }

    /* find indexes belonging to the last index set */
    find_indexes(index_sets, n1*n2, path[n1*n2-1], indexes);
    for (j = 0; j < num_indexes; j++)
        index_sets_new[indexes[j]-1] = floor((n-1)/2)+1;

    return path;

}

float **reorder(float **data,   /* data array */
               int *path,       /* path vector */
               int nt,          /* length of time samples */
               int n1, int n2,  /* size of the index sets */
               bool inv         /* inverse transform */)
/*< reorder the data according to the path vector >*/
{
    int i;
    float **data_new     = sf_floatalloc2(nt, n1*n2);   // new data array corresponding to the path 

    if (inv) {
        for (i = 0; i < n1*n2; i++)
            for (j = 0; j < nt; j++)
                data_new[path[i]-1][j] = data[i][j];
    } else {
        for (i = 0; i < n1*n2; i++) 
            for (j = 0; j < nt; j++)
                data_new[i][j] = data[path[i]-1][j];
    }

    return data_new;
}
