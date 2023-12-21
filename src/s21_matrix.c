#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  if (rows <= 0 || columns <= 0 || result == NULL) return INVALID_MATRIX;
  int ext = OK;
  result->rows = rows;
  result->columns = columns;
  result->matrix = calloc(rows, sizeof(double *));
  if (result->matrix == NULL) {
    ext = CALCULATION_ERROR;
  } else {
    for (int i = 0; i < rows; i++) {
      result->matrix[i] = calloc(columns, sizeof(double));
      if (result->matrix[i] == NULL) {
        s21_remove_matrix(result);
        ext = CALCULATION_ERROR;
        break;
      }
    }
  }
  return ext;
}

void s21_remove_matrix(matrix_t *A) {
  if (A->matrix) {
    for (int i = 0; i < A->rows; i++) {
      free(*(A->matrix + i));
      *(A->matrix + i) = NULL;
    }
    free(A->matrix);
  }
  A->rows = 0;
  A->columns = 0;
  A->matrix = NULL;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  if (A->rows != B->rows || A->columns != B->columns) return FAILURE;
  int ext = SUCCESS;
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      double diff = fabs(A->matrix[i][j] - B->matrix[i][j]);
      if (diff > EPSILON) {
        ext = FAILURE;
        break;
      }
    }
  }
  return ext;
}

void s21_sum(matrix_t *A, matrix_t *B, matrix_t *result) {
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
    }
  }
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL) return INVALID_MATRIX;
  int ext = OK;
  if (A->rows != B->rows || A->columns != B->columns) {
    ext = CALCULATION_ERROR;
  } else {
    ext = s21_create_matrix(A->rows, A->columns, result);
    if (ext == OK) {
      s21_sum(A, B, result);
    }
  }
  return ext;
}

void s21_sub(matrix_t *A, matrix_t *B, matrix_t *result) {
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
    }
  }
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL) return INVALID_MATRIX;
  int ext = OK;
  if (A->rows != B->rows || A->columns != B->columns) {
    ext = CALCULATION_ERROR;
  } else {
    ext = s21_create_matrix(A->rows, A->columns, result);
    if (ext == OK) {
      s21_sub(A, B, result);
    }
  }
  return ext;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  if (A == NULL || result == NULL) return INVALID_MATRIX;
  int ext = s21_create_matrix(A->rows, A->columns, result);
  if (ext == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  }
  return ext;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  if (A == NULL || B == NULL || result == NULL) return INVALID_MATRIX;
  int ext = OK;
  if (A->rows == B->columns) {
    ext = s21_create_matrix(A->rows, B->columns, result);
    if (ext == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < B->columns; j++) {
          result->matrix[i][j] = 0;
          for (int k = 0; k < A->columns; k++) {
            result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
          }
        }
      }
    }
  } else {
    ext = CALCULATION_ERROR;
  }
  return ext;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  if (A == NULL || result == NULL) return INVALID_MATRIX;
  int ext = s21_create_matrix(A->columns, A->rows, result);
  if (ext == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[j][i] = A->matrix[i][j];
      }
    }
  }
  return ext;
}

void get_minor(int row, int col, matrix_t *A, matrix_t *result) {
  int act_row = 0, act_col = 0;
  for (int i = 0; i < A->rows; i++) {
    if (i == row) continue;
    act_col = 0;
    for (int j = 0; j < A->rows; j++) {
      if (j == col) continue;
      result->matrix[act_row][act_col] = A->matrix[i][j];
      act_col++;
    }
    act_row++;
  }
}

double get_determinant(matrix_t *A) {
  double det = 0;
  if (A->rows == 1) {
    det = A->matrix[0][0];
  } else if (A->rows == 2) {
    det = A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
  } else {
    double sign = 1;
    matrix_t minor;
    int ext = OK;
    for (int i = 0; i < A->rows && ext == OK; i++) {
      ext = s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
      if (ext == OK) {
        get_minor(0, i, A, &minor);
        det += A->matrix[0][i] * sign * get_determinant(&minor);
        sign *= -1;
        s21_remove_matrix(&minor);
      }
    }
  }
  return det;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  if (A == NULL || result == NULL || A->rows <= 0 || A->columns <= 0)
    return INVALID_MATRIX;
  int ext = OK;
  if (A->rows == A->columns && A->rows == 1) {
    s21_create_matrix(1, 1, result);
    result->matrix[0][0] = 1;
  } else if (A->rows == A->columns && A->rows > 1) {
    s21_create_matrix(A->rows, A->columns, result);
    for (int i = 0; i < A->rows && ext == OK; i++) {
      for (int j = 0; j < A->columns && ext == OK; j++) {
        matrix_t minor;
        ext = s21_create_matrix(A->rows - 1, A->columns - 1, &minor);
        if (ext == OK) {
          get_minor(i, j, A, &minor);
          result->matrix[i][j] = pow(-1, i + j) * get_determinant(&minor);
        }
        s21_remove_matrix(&minor);
      }
    }
  } else {
    ext = CALCULATION_ERROR;
  }
  return ext;
}

int s21_determinant(matrix_t *A, double *result) {
  if (A == NULL) return INVALID_MATRIX;
  if (A->rows != A->columns) return CALCULATION_ERROR;
  int ext = OK;
  if (A->rows > 0) {
    *result = get_determinant(A);
  } else {
    ext = INVALID_MATRIX;
  }
  return ext;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  double det = 0;
  int ext = s21_determinant(A, &det);
  if (ext == OK && det != 0) {
    matrix_t calc_com = {0};
    ext = s21_calc_complements(A, &calc_com);
    if (ext == OK) {
      matrix_t transpose = {0};
      ext = s21_transpose(&calc_com, &transpose);
      s21_remove_matrix(&calc_com);
      if (ext == OK) {
        s21_mult_number(&transpose, 1 / det, result);
        s21_remove_matrix(&transpose);
      }
    }
  }
  return det == 0 && ext == 0 ? CALCULATION_ERROR : ext;
}
