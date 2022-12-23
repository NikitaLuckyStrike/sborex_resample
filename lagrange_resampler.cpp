//
// Created by n.titov on 17.11.2022.
//

#include "lagrange_resampler.h"


lagrange_resampler::lagrange_resampler(float t0, int size_old, float *t_old, float *v_old, float dt_new)
    //: t0(t0), size_old(size_old), t_old(t_old), v_old(v_old), dt_new(dt_new)
{
  set_data(t_old, v_old);
  set_param(t0, size_old, dt_new);

  if (t_old == nullptr || v_old == nullptr || t0 <=0 || size_old <= 0 || dt_new <= 0) {
    this->has_data = SOME_ERROR;
  } else {
    this->new_size = calc_new_size(this->dt_new, this->size_old, this->t_old);
    this->t_new = new float[new_size]; //трассы после ресэмплинга
    this->v_new = new float[new_size];
    this->has_data = GOT_DATA;
  } //add set_data() here
}

lagrange_resampler::lagrange_resampler() { //конструктор без инициализации переменных
}



void lagrange_resampler::set_param(float t0, int size_old, float dt_new) {
  if (t0 <= 0 || size_old <= 0 || dt_new <= 0) {
    this->has_data = SOME_ERROR;
  } else {
    this->t0 = t0;
    this->size_old = size_old;
    this->dt_new = dt_new;
    this->new_size = calc_new_size(this->dt_new, this->size_old, this->t_old);
  }
}

void lagrange_resampler::set_data(float *t_old, float *v_old) {
  if (t_old == nullptr || v_old == nullptr) {
    this->has_data = SOME_ERROR;
  } else {
    this->t_old = t_old;
    this->v_old = v_old;
    this->has_data = GOT_DATA;
  }
}

// функция, позволяющая передать объекту значения, необходимые для работы
//если в объекте уже есть данные, то они перезапишутся
void lagrange_resampler::set_data(float t0, int size_old, float* t_old, float* v_old, float dt_new) {
  if (t_old == nullptr || v_old == nullptr || t0 <=0 || size_old <= 0 || dt_new <= 0) {
    this->has_data = SOME_ERROR;
  } else {
    if (this->has_data == GOT_DATA || this->has_data == READY_DATA) {
      clear();
    }
    
    this->has_data = GOT_DATA;
    
    this->t0 = t0;
    this->size_old = size_old;
    if (t_old == nullptr) {
      this->has_data = SOME_ERROR;
    } else {
      this->t_old = t_old;
    }
      if (v_old == nullptr) {
      this->has_data = SOME_ERROR;
    } else {
      this->v_old = v_old;
    }

    
    this->new_size = calc_new_size(this->dt_new, this->size_old, this->t_old);
    this->t_new = new float[new_size]; //трассы после ресэмплинга
    this->v_new = new float[new_size];
  }
}

lagrange_resampler::~lagrange_resampler(){
  clear();
}

int lagrange_resampler::get_new_size() {
  return this->new_size;
}

int lagrange_resampler::calc_new_size(float dt_new, int old_size, float* old_t) {
  return ((old_t[old_size-1] - old_t[0])/ dt_new + ((old_t[old_size-1] - old_t[0])/old_size)/ dt_new);
}

float lagrange_resampler::lagrange_point(float x, int n, float *t_arr, float *y_arr) {
  // x - new time;
  // n - lagrange polinominal size (передаем пока 4)
  // t_arr - array of time values
  // y_arr - array of values
  float sum = 0;
  for (int i = 0; i < n; ++i) {
    float l = 1; // многочлен Лагранжа
    for (int j = 0; j < n; ++j)
      if (j != i)
        l *= (x - t_arr[j]) / (t_arr[i] - t_arr[j]);
    sum += y_arr[i] * l;
  }
  return sum;
}

int lagrange_resampler::find_index_for_lagrange(float *t_old, int size_old, float t_cur,
                            int old_t_cur_ind) {
  int result = -1;
  for (int i = old_t_cur_ind; i < size_old; i++) { //рассчитать рамки
    if (t_cur >= t_old[i + 1] && t_cur <= t_old[i + 2]) {
      result = i;
      break;
    }
  }
  return result;
}

int lagrange_resampler::lagrange_trace(float t0, float dt_new, int size_old, float *t_old, float *v_old,
                   int size_new, float *t_new, float *v_new) {

  // t0 - время синхроимпульса, dt -новый шаг дискретизации, size - размер
  // старого массива, t_old v_old - старая трасса, t_new v_new - новая трасса

  //int size_new = 0; //размер трассы после ресамплинга
  float t = t0; //переменная для рассчета длины новой трассы
  int old_t_cur_ind = 0; //переменная для хранения индекса массива отсчетов
  //времени старой трассы

  int ind_t_new = 0; // индекс для элементов новой трассы для этого цикла
//  while (t <= t_old[size_old - 1]) { // узнаем кол-во отсчетов в новой трассе
//    // и заполняем элементами массив временных отсчетов для новой трассы
//    size_new++;
//    t_new[ind_t_new] = t;
//    ind_t_new++;
//    t += dt;
//  }

  while (ind_t_new<size_new) {
    t_new[ind_t_new] = t;
    t+=dt_new;
    ind_t_new++;
  }

  for (int ind_t_new = 0; ind_t_new < size_new; ind_t_new++) {

    if (t_new[ind_t_new] < t_old[1]) { //для краевого случая, когда отсчет новой
      //трассы до 2-го отсчета старой

      float tmp_t_old[4]; //болванки для расчета
      float tmp_v_old[4];
      for (int j = 0; j < 4; j++) {
        tmp_t_old[j] = t_old[j];
        tmp_v_old[j] = v_old[j];
      }
      v_new[ind_t_new] = lagrange_point(t_new[ind_t_new], LAGRANGE_DEGREE,
                                        tmp_t_old, tmp_v_old);

    } else if (t_new[ind_t_new] > t_old[size_old - 2]) {
      //для краевого случая, когда отсчет новой трассы
      //позже предпоследнего отсчета старой

      float tmp_t_old[4]; //болванки для расчета
      float tmp_v_old[4];
      for (int j = 0; j < 4; j++) {
        tmp_t_old[j] = t_old[size_old - 4 + j];
        tmp_v_old[j] = v_old[size_old - 4 + j];
      }
      v_new[ind_t_new] = lagrange_point(t_new[ind_t_new], LAGRANGE_DEGREE,
                                        tmp_t_old, tmp_v_old);
    } else {
      // для всех остальных случаев

      old_t_cur_ind = find_index_for_lagrange(t_old, size_old, t_new[ind_t_new],
                                              old_t_cur_ind);

      float tmp_t_old[4]; //болванки для расчета
      float tmp_v_old[4];
      for (int j = 0; j < 4; j++) {
        tmp_t_old[j] = t_old[old_t_cur_ind + j];
        tmp_v_old[j] = v_old[old_t_cur_ind + j];
      }
      v_new[ind_t_new] = lagrange_point(t_new[ind_t_new], LAGRANGE_DEGREE,
                                        tmp_t_old, tmp_v_old);
    }

  }
  return ind_t_new;
}

void lagrange_resampler::clear() {
  this->has_data = 0;

  this->dt_new = 0;

  this->size_old = 0; //LEN
  this->new_size = 0;

  this->t_old = nullptr;
  this->v_old = nullptr;

  delete t_new;
  delete v_new;


  this->t0 = 0;
}

int lagrange_resampler::do_approximation() {
  if (has_data == GOT_DATA) {
    lagrange_trace(this->t0, this->dt_new, this->size_old, this->t_old, this->v_old, this->new_size,
                 this->t_new, this->v_new);
  }
  has_data = READY_DATA;
  return has_data;
}
float* lagrange_resampler::get_new_v() {
  if (has_data == READY_DATA) {
    return this->v_new;
  } else {
    return nullptr;
  }
}
float* lagrange_resampler::get_new_t() {
  if (has_data == READY_DATA) {
    return this->t_new;
  } else {
    return nullptr;
  }
}