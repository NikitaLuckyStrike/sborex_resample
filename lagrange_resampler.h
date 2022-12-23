//
// Created by n.titov on 17.11.2022.
//

#ifndef LAGRANGE_CPP_LAGRANGE_RESAMPLER_H
#define LAGRANGE_CPP_LAGRANGE_RESAMPLER_H

#include "resampler.h"


#define LAGRANGE_DEGREE 4

#define NO_DATA 0
#define GOT_DATA 1
#define READY_DATA 2
#define SOME_ERROR -1


class lagrange_resampler : public resampler {

public:
  lagrange_resampler(float t0, int size_old, float* t_old, float* v_old, float dt_new);
  lagrange_resampler();
  ~lagrange_resampler();

  int do_approximation(); // override; //апроксимирует значения, пишет новые значения в массивы t_new v_new
  void clear(); // очистка всех имеющихся
  float* get_new_v(); //геттеры для новой трассы
  float* get_new_t();

  void set_data(float t0, int size_old, float* t_old, float* v_old, float dt_new);
  void set_data(float* t_old, float* v_old);
  void set_param(float t0, int size_old, float dt_new);

  int get_new_size();

private:
  int lagrange_trace( float t0, float dt, int size_old, float *t_old,
                     float *v_old,int size_new, float *t_new, float *v_new);
  //функция для рассчета значений трасы после ресамплинга

  float* v_new; // arrays for trace
  float* t_new;

  float dt_new;

  int size_old; //LEN
  int new_size;

  float* v_old = nullptr;
  float* t_old = nullptr;


  int has_data = NO_DATA;


  float t0;

  float lagrange_point(float x, int n, float *t_arr,
      float *y_arr); // рассчет значений для каждого отдельного дискрета

  int find_index_for_lagrange(float *t_old, int size_old, float t_cur,
  int old_t_cur_ind);// поиск индексов элементов для рассчета многочлена Лагранжа
// для каждого отсчета времени новой трассы
  int calc_new_size(float dt_new, int old_size, float* old_t);
};

#endif // LAGRANGE_CPP_LAGRANGE_RESAMPLER_H
