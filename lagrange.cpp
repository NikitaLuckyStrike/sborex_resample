
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <thread> //потоки
#include <chrono> //время

#include "lagrange_resampler.h"

#define LEN 8000
#define LAGRANGE_DEGREE 4  // СТЕПЕНЬ МНОГОЧЛЕНА ЛАГРАНЖА

float lagrange_point(
    float x, int n, float *t_arr,
    float *y_arr);  // рассчет значений для каждого отдельного дискрета

int lagrange_trace(
    float t0, float dt, int size_old, float *t_old, float *v_old, int size_new,
    float *t_new,
    float *v_new);  //функция для рассчета значений трасы после ресамплинга

void generate_trace(
    float *time, float *value, float start_time, float dt,
    int trace_lenght);  // создание трассы с синусом для отладки алгоритма

void write_output(float *t, float *v, int size,
                  char *filename);  //вывод в текстовый файл для проверки

int find_index_for_lagrange(float *t_old, int size_old, float t_cur,
                            int old_t_cur_ind);
// поиск индексов элементов для рассчета многочлена Лагранжа
// для каждого отсчета времени новой трассы

int calc_new_size(float new_dt, int old_size, float *old_t);

using namespace std;

//int main(int argc, char** argv) {
//#ifndef LEN
//#define LEN 8000
//#endif
//  float t[LEN] = {0};  // параметры трассы для отладки
//  float v[LEN] = {0};
//  float start_time = 1234.5;
//  float old_dt = 0.5;
//  float new_dt = 0.1;
//
//  generate_trace(t, v, start_time, old_dt, LEN);  // test trace
//  lagrange_resampler *test =
//      new lagrange_resampler(start_time, 8000, t, v, new_dt);
//  write_output(t, v, LEN,
//               "C:\\Users\\Lenovo\\Desktop\\sborex_resample\\out_old.txt");
//  test->do_approximation();
//  float *tf = test->get_new_t();
//  float *vf = test->get_new_v();
//  write_output(test->get_new_t(), test->get_new_v(), test->get_new_size(),
//               "C:\\Users\\Lenovo\\Desktop\\sborex_resample\\out_new.txt");
//  test->clear();
//  return 0;
//}

#define NOT 16 //number of traces
int main() {

  float t_big[LEN*NOT] = {0};
  float v_big[LEN*NOT] = {0};
  float start_time = 1234.5;
  float old_dt = 0.5;
  float new_dt = 0.1;

  for (int i = 0; i < NOT; i++) {
      generate_trace(t_big+LEN*i,v_big+LEN*i,start_time, old_dt, LEN);
  }
  printf("traces generated\n");
    return 0;
}









int calc_new_size(float new_dt, int old_size, float *old_t) {
  return ((old_t[old_size - 1] - old_t[0]) / new_dt +
          ((old_t[old_size - 1] - old_t[0]) / old_size) / new_dt);
}

float lagrange_point(float x, int n, float *t_arr, float *y_arr) {
  // x - new time;
  // n - lagrange polinominal size (передаем пока 4)
  // t_arr - array of time values
  // y_arr - array of values
  float sum = 0;
  for (int i = 0; i < n; ++i) {
    float l = 1;  // многочлен Лагранжа
    for (int j = 0; j < n; ++j)
      if (j != i) l *= (x - t_arr[j]) / (t_arr[i] - t_arr[j]);
    sum += y_arr[i] * l;
  }
  return sum;
}

int lagrange_trace(float t0, float dt_new, int size_old, float *t_old,
                   float *v_old, int size_new, float *t_new, float *v_new) {
  // t0 - время синхроимпульса, dt -новый шаг дискретизации, size - размер
  // старого массива, t_old v_old - старая трасса, t_new v_new - новая трасса

  // int size_new = 0; //размер трассы после ресамплинга
  float t = t0;  //переменная для рассчета длины новой трассы
  int old_t_cur_ind = 0;  //переменная для хранения индекса массива отсчетов
                          //времени старой трассы

  int ind_t_new = 0;  // индекс для элементов новой трассы для этого цикла
  //  while (t <= t_old[size_old - 1]) { // узнаем кол-во отсчетов в новой
  //  трассе
  //    // и заполняем элементами массив временных отсчетов для новой трассы
  //    size_new++;
  //    t_new[ind_t_new] = t;
  //    ind_t_new++;
  //    t += dt;
  //  }

  while (ind_t_new < size_new) {
    t_new[ind_t_new] = t;
    t += dt_new;
    ind_t_new++;
  }

  for (int ind_t_new = 0; ind_t_new < size_new; ind_t_new++) {
    if (t_new[ind_t_new] < t_old[1]) {  //для краевого случая, когда отсчет
                                        //новой трассы до 2-го отсчета старой

      float tmp_t_old[4];  //болванки для расчета
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

      float tmp_t_old[4];  //болванки для расчета
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

      float tmp_t_old[4];  //болванки для расчета
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

int find_index_for_lagrange(float *t_old, int size_old, float t_cur,
                            int old_t_cur_ind) {
  int result = -1;
  for (int i = old_t_cur_ind; i < size_old; i++) {  //рассчитать рамки
    if (t_cur >= t_old[i + 1] && t_cur <= t_old[i + 2]) {
      result = i;
      break;
    }
  }
  return result;
}

void generate_trace(float *time, float *value, float start_time, float dt,
                    int trace_lenght) {
  //    std::srand(std::time(nullptr));
  for (int i = 0; i < LEN; i++) {
    time[i] = start_time + dt * i;
    value[i] = sin(time[i] / 20) * 100;
    // printf("%d \t time = %.2lf\tvalue = %.2lf\n", i, (time[i]), (value[i]));
  }
}

void write_output(float *t, float *v, int size, char *filename) {
  ofstream file(filename);
  //"C:\\Users\\n.titov\\Desktop\\sborex_resample\\out.txt"
  //  float dt = trace_size_new>trace_size_old? t_new[1] - t_new[0]:t_old[1] -
  //  t_old[0]; // dt of new
  for (int i = 0; i < size; i++) {
    file << t[i] << "\t" << v[i] << "\n";
  }
  file.close();
}


