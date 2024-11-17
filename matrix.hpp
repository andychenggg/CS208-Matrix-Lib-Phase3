#ifndef CPPPROJECTS_MATRIX_HPP
#define CPPPROJECTS_MATRIX_HPP

#include <iostream>
#include <memory.h>
#include <vector>
#include "memory"
#include "initializer_list"
#include "MATexception.hpp"
#include "random"
#include <minmax.h>

#define __WHERE__ __FILE__, __LINE__, __func__

//编译优化
#pragma GCC optimize(3, "Ofast", "inline")
#pragma GCC optimize(1)
#pragma GCC optimize(2)
#pragma GCC optimize(3)
//使用到的std命名空间的一些变量
using std::cout;
using std::cerr;
using std::endl;
using std::vector;

namespace mat {
    //支持的特殊矩阵的类型
    enum mat_type {
        identity_mat = 0,   //单位矩阵
    };
    //支持两种打印模式
    enum io_mode {
        channel_by_channel = 0,   //每一次将一个通道全部进行输入或者输出
        row_by_row_col_by_col = 1 //每一次将一个位置的所有通道进行输入或者输出
    };

    template<typename T>
    class matrix;

    template<typename T>
    class mat_iterator;

    class mat_range {
    private:
        //上下界
        long long lb, ub;  //表示有效范围在[lb, ub)
    public:
        template<class T>
        friend
        class matrix;

        mat_range() {}

        mat_range(long long lb, long long ub) : lb(lb), ub(ub) {}

        mat_range(const mat_range &m1) {
            ub = m1.ub;
            lb = m1.lb;
        }

        operator long long() const {
            return ub - lb;
        }

        [[nodiscard]] long long getUB() const {
            return ub;
        }


        [[nodiscard]] long long getLB() const {
            return lb;
        }

    };

    template<typename T>
    class matrix {
    private:
        size_t rows;
        size_t cols;
        vector <std::shared_ptr<T>> data;
        size_t channels;
        //为了兼容ROI的特性
        size_t shift;   //实际数据位置距离当层0号元素的偏移量
        size_t step;    //父矩阵的列数
        size_t step_row;//父矩阵的行数
        mat_range channel_range;
        //迭代器
        mat_iterator<T> begin;  //头指针
        mat_iterator<T> cur;    //当前指针
        mat_iterator<T> end;    //尾指针的下一个指针
        //输入输出方式
        io_mode iom = io_mode::channel_by_channel;

        //自定义数据析构函数
        static void self_def_deleter(T *pt) {
            delete[] pt;
        }

    public:
        //友元类
        friend class mat_iterator<T>;

        template<class V>
        friend
        class matrix;
        /*
         * 1、构造器与析构函数
         */
        //不允许隐式转换
        explicit matrix(size_t rows = 1, size_t cols = 1, size_t channels = 1) noexcept(true);

        //使用已有数据进行初始化
        matrix(size_t rows, size_t cols, size_t channels, const T *data) noexcept(true);

        //使用初始化列表进行初始化
        matrix(size_t rows, size_t cols, size_t channels, const std::initializer_list<T> &list) noexcept(true);

        //全部设为一个值
        matrix(size_t rows, size_t cols, size_t channels, T t) noexcept(true);

        //创建特殊矩阵
        matrix(size_t rows, size_t channels, mat_type type) noexcept(true);

        //拷贝构造器
        matrix(const matrix &m) noexcept(true);

        //析构函数
        ~matrix();

        /*
         * 1、ROI
         */
        //ROI截取构造器
        matrix(matrix<T> &m, mat_range rowRange, mat_range colRange, mat_range chanRange) noexcept(false);

        //增加一层矩阵
        void add_channel(std::initializer_list<T> list) noexcept;

        //删除指定层的矩阵
        void erase_channel(size_t pos) noexcept;

        //ROI矩阵变大变小
        void reshape(long long row, long long col, long long chan) noexcept(false);

        //ROI矩阵上下左右前后的移动
        void move(long long row, long long col, long long chan) noexcept(false);

        //判断是否是子矩阵
        bool isSubMatrix(const matrix<T> &m1) const noexcept;

        /*
         * 3、矩阵的setter和getter方法：部分使用运算符重载
         */
        //赋值操作
        matrix<T> &operator=(const matrix<T> &matrix1) noexcept;

        //硬拷贝
        void hardCopy(const matrix<T> &m1) noexcept;

        //访问第chan层的矩阵头指针，只能是右值
        T *operator[](size_t chan) const noexcept(false);

        //访问第i层j行k列的元素,可赋值
        T &operator()(size_t i, size_t j, size_t k) const noexcept(false);

        //返回行数
        [[nodiscard]] size_t getRows() const noexcept;

        //返回列数
        [[nodiscard]] size_t getCols() const noexcept;

        //返回通道数
        [[nodiscard]] size_t getChannels() const noexcept;

        //返回父矩阵的列数
        [[nodiscard]] size_t getStep() const noexcept;

        //返回父矩阵的行数
        [[nodiscard]] size_t getStepRows() const noexcept;

        //返回偏移量
        [[nodiscard]] size_t getShift() const noexcept;


        //setter
        //将第chan个通道全部设为一个值
        void set_channel(size_t chan, T val) noexcept(false);

        //将第chan个通道第i行设为一个值
        void set_row(size_t chan, size_t i, T val) noexcept(false);

        //将第chan个通道第j列设为一个值
        void set_col(size_t chan, size_t j, T val) noexcept(false);

        //迭代器
        mat_iterator<T> getBegin() {
            return begin;
        }

        mat_iterator<T> getEnd() {
            return end;
        }


        /*
         * 4、运算符重载
         */
        /*
         * 加法
         */
        //矩阵和常数相加
        template<typename V>
        [[nodiscard]]  //开辟了新内存，不应舍弃
        matrix<T> operator+(V val) const noexcept;

        template<typename V>
        [[nodiscard]]  //开辟了新内存，不应舍弃
        friend matrix<T> operator+(V val, const matrix<T> &m) noexcept {
            return m + val;
        }

        //矩阵和矩阵相加
        template<class V>
        [[nodiscard]]  //开辟了新内存，不应舍弃
        matrix<T> operator+(const matrix<V> &m) const noexcept(false);

        //运算符+=
        template<typename V>
        matrix<T> &operator+=(V val) noexcept;

        matrix<T> &operator+=(const matrix<T> &m) noexcept(false);

        /*
         * 乘法
         */
        //矩阵和常数相乘
        template<typename V>
        [[nodiscard]]  //开辟了新内存，不应舍弃
        matrix<T> operator*(V val) const noexcept;

        template<typename V>
        [[nodiscard]]  //开辟了新内存，不应舍弃
        friend matrix<T> operator*(V val, const matrix<T> &m) noexcept {
            return m * val;
        }

        //矩阵和矩阵相乘
        [[nodiscard]]  //开辟了新内存，不应舍弃
        matrix<T> operator*(const matrix<T> &m) const noexcept(false);

        /*
         * 减法
         */
        //矩阵和常数相减
        template<typename V>
        [[nodiscard]]  //开辟了新内存，不应舍弃
        matrix<T> operator-(V val) const noexcept;

        template<typename V>
        [[nodiscard]]  //开辟了新内存，不应舍弃
        friend matrix<T> operator-(V val, const matrix<T> &m) noexcept {
            return (-1) * m + val;
        }

        //矩阵和矩阵相减
        [[nodiscard]]  //开辟了新内存，不应舍弃
        matrix<T> operator-(const matrix<T> &m) const noexcept(false);

        //运算符-=
        template<typename V>
        matrix<T> &operator-=(V val) noexcept;

        matrix<T> &operator-=(const matrix<T> &m) noexcept(false);

        /*
         * 按位运算符：在这里只重载 ^ & | ~ \<\< >> 这几个最基本的按位运算符的最基本形式
         */
        //矩阵与数之间的^ & |
        template<typename V>
        [[nodiscard]]  //开辟了新内存，不应舍弃
        matrix<T> operator^(V val) const noexcept;

        template<typename V>
        [[nodiscard]]  //开辟了新内存，不应舍弃
        matrix<T> operator&(V val) const noexcept;

        template<typename V>
        [[nodiscard]]  //开辟了新内存，不应舍弃
        matrix<T> operator|(V val) const noexcept;

        //矩阵与矩阵之间的按位^ & |
        [[nodiscard]]  //开辟了新内存，不应舍弃
        matrix<T> operator^(const matrix<T> &m) const noexcept(false);

        [[nodiscard]]  //开辟了新内存，不应舍弃
        matrix<T> operator&(const matrix<T> &m) const noexcept(false);

        [[nodiscard]]  //开辟了新内存，不应舍弃
        matrix<T> operator|(const matrix<T> &m) const noexcept(false);




        /*
         * 比较运算符的重载
         */
        //==: 矩阵规模和每个矩阵元素必须相等
        bool operator==(const matrix<T> &m) const noexcept(false);

        //\!=: 规模不同或至少有一个元素不相等
        bool operator!=(const matrix<T> &m) const noexcept(false);

        // >=: 规模相同的情况下，每个元素都要大于等于， >: 每个元素都要大于
        bool operator>=(const matrix<T> &m) const noexcept(false);

        bool operator>(const matrix<T> &m) const noexcept(false);

        // <=: 每个元素都要大于等于， <: 每个元素都要大于
        bool operator<=(const matrix<T> &m) const noexcept(false);

        bool operator<(const matrix<T> &m) const noexcept(false);

        /*
         * 自增, 自减运算符；
         */
        //++，表示每个位置自增1；
        matrix<T> &operator++() noexcept {
            (*this) += 1;
            return *this;
        }

        matrix<T> operator++(int) noexcept {
            matrix<T> m1(rows, cols, channels);
            m1.hardCopy(*this);
            operator++();
            return m1;
        }

        //矩阵自身移位操作
        //右移
        template<typename V>
        matrix<T> &operator<<=(V val);

        template<typename V>
        matrix<T> operator<<(V val);

        //左移
        template<typename V>
        matrix<T> &operator>>=(V val);

        template<typename V>
        matrix<T> operator>>(V val);


        /*
         * 输入，输出运算符
         */
        //设置io模式
        void set_IO_mode(io_mode m) {
            iom = m;
        }
        //输出

        friend
        std::ostream &operator<<(std::ostream &o1, matrix<T>& m1) {
            if (m1.iom == io_mode::channel_by_channel) {
                for (int i = 0; i < m1.rows; i++) {
                    if (i == 0) {
                        o1 << "[ ";
                    } else {
                        o1 << "  ";
                    }
                    for (int j = 0; j < m1.cols; j++) {
                        for (size_t k = m1.channel_range.getLB(); k < m1.channel_range.getUB(); k++) {
                            if (k == m1.channel_range.getLB()) {
                                o1 << "{";
                            }
                            o1 << m1.data[k].get()[m1.shift + i * m1.step + j];
                            if (k == m1.channel_range.getUB() - 1) {
                                o1 << "}";
                            }
                            if (i != m1.rows - 1 || j != m1.cols - 1 || k != m1.channel_range.getUB() - 1) {
                                o1 << ", ";
                            }
                        }
                    }
                    if (i != m1.rows - 1) {
                        o1 << endl;
                    } else {
                        o1 << " ]" << endl;
                    }
                }
            } else {
                for (size_t k = m1.channel_range.getLB(); k < m1.channel_range.getUB(); k++) {
                    if (k == m1.channel_range.getLB()) {
                        o1 << "[";
                    } else {
                        o1 << " ";
                    }
                    for (int i = 0; i < m1.rows; i++) {
                        if (i == 0) {
                            o1 << "{";
                        }
                        else
                        {
                            o1 << "  ";
                        }
                        for (int j = 0; j < m1.cols; j++) {
                            o1 << m1.data[k].get()[m1.shift + i * m1.step + j];
                            if (i != m1.rows - 1 || j != m1.cols - 1) {
                                o1 << ", ";
                            }
                        }
                        if (i == m1.rows - 1) {
                            o1 << "}";
                        }
                        else{
                            cout<<endl;
                        }
                    }
                    if (k != m1.channel_range.getUB() - 1) {
                        o1 << endl;
                    } else {
                        o1 << "]" << endl;
                    }
                }
            }
            return o1;
        }

        //输入
        friend
        std::istream &operator>>(std::istream &i1, matrix<T>& m1) {
            if (m1.iom == io_mode::channel_by_channel) {
                for (int i = 0; i < m1.rows; i++) {
                    for (int j = 0; j < m1.cols; j++) {
                        for (size_t k = m1.channel_range.getLB(); k < m1.channel_range.getUB(); k++) {
                            i1 >> m1.data[k].get()[m1.shift + i * m1.step + j];
                            }
                    }
                }
            } else {
                for (size_t k = m1.channel_range.getLB(); k < m1.channel_range.getUB(); k++) {

                    for (int i = 0; i < m1.rows; i++) {
                        for (int j = 0; j < m1.cols; j++) {
                            i1 >>  m1.data[k].get()[m1.shift + i * m1.step + j];
                        }
                    }

                }
            }
            return i1;
        }
    };


    template<typename T>
    class mat_iterator {
    private:
        //考虑到ROI
        //父矩阵的数据vector
        vector <std::shared_ptr<T>> begin;
        //当前所在的位置
        T *current;
        //可以移动的绝对范围
        mat_range row_range, col_range, channel_range;
        //目前所在的绝对位置
        size_t cur_row, cur_col, cur_chan;
        size_t shift, step, step_row;
    public:
        friend class matrix<T>;

        /*
         * 构造器
         */
        mat_iterator() {}

        mat_iterator(const matrix<T> &m1, T *current, size_t i_in_m1, size_t j_in_m1, size_t k_in_m1) {
            begin = m1.data;
            this->current = current;
            row_range = mat_range(m1.shift % m1.step, m1.shift % m1.step + m1.rows);
            col_range = mat_range(m1.shift / m1.step, m1.shift % m1.step + m1.cols);
            channel_range = m1.channel_range;
            cur_row = (row_range.getLB() + i_in_m1);
            cur_col = (col_range.getLB() + j_in_m1);
            cur_chan = channel_range.getLB() + k_in_m1;
            shift = m1.shift;
            step = m1.step;
            step_row = m1.step_row;
        }

        mat_iterator(const mat_iterator &m1) :
                begin(m1.begin),
                current(m1.current),
                row_range(m1.row_range), col_range(m1.col_range), channel_range(m1.channel_range),
                cur_row(m1.cur_row), cur_col(m1.cur_col), cur_chan(m1.cur_chan),
                shift(m1.shift), step(m1.step), step_row(m1.step_row) {}

        /*
         * 自增自减运算符
         */
        inline mat_iterator<T> &operator++() noexcept(false) {
            if (cur_row == row_range.getUB() - 1 &&
                cur_col == col_range.getUB() &&
                cur_chan == channel_range.getUB() - 1) {
                IteratorOutOfBoundException e("Current is equal or bigger than end!", __WHERE__);
#ifdef __NOEXCEPT_TRUE__
                cout << endl;
                cerr << e.what();
                return *this;
#else
                throw e;
#endif
            }
                //可以到后一格，即与end相同
            else if (cur_row == row_range.getUB() - 1 &&
                     cur_col == col_range.getUB() - 1 &&
                     cur_chan == channel_range.getUB() - 1) {
                cur_col++;
                current++;
                return *this;
            } else if (cur_row == row_range.getUB() - 1 &&
                       cur_col == col_range.getUB() - 1) {
                cur_chan++;
                cur_row = row_range.getLB();
                cur_col = col_range.getLB();
                current = &(begin[cur_chan].get()[shift]);
                return *this;
            } else if (cur_col == row_range.getUB() - 1) {
                cur_row++;
                cur_col = col_range.getLB();
                current = current + 1 + step - col_range;
                return *this;
            } else {
                cur_col++;
                current++;
                return *this;
            }
        }

        inline mat_iterator<T> operator++(int) noexcept(false) {
            mat_iterator<T> it = *this;
            operator++();
            return it;
        }

        inline mat_iterator<T> &operator--() noexcept(false) {
            if (cur_row == row_range.getLB() &&
                cur_col == col_range.getLB() &&
                cur_chan == channel_range.getLB()) {
                IteratorOutOfBoundException e("Current is equal or bigger than end!", __WHERE__);
#ifdef __NOEXCEPT_TRUE__
                cout << endl;
                cerr << e.what();
                return *this;
#else
                throw e;
#endif
            } else if (cur_row == row_range.getLB() &&
                       cur_col == col_range.getLB()) {
                cur_chan--;
                cur_row = row_range.getUB() - 1;
                cur_col = col_range.getUB() - 1;
                current = &(begin[cur_chan].get()[shift + (row_range - 1) * step + col_range - 1]);
                return *this;
            } else if (cur_col == row_range.getLB()) {
                cur_row--;
                cur_col = col_range.getUB() - 1;
                current = current - 1 - step - col_range;
                return *this;
            } else {
                cur_col--;
                current--;
                return *this;
            }
        }

        inline mat_iterator<T> operator--(int) noexcept(false) {
            mat_iterator<T> it = *this;
            operator--();
            return it;
        }

        /*
         * 赋值运算符
         */
        mat_iterator<T> &operator=(const mat_iterator<T> &it) noexcept {
            if (&it == this) {
                return *this;
            } else {
                begin = it.begin;
                current = it.current;

                row_range = it.row_range;
                col_range = it.col_range;
                channel_range = it.channel_range;

                cur_row = it.cur_row;
                cur_col = it.cur_col;
                cur_chan = it.cur_chan;

                shift = it.shift;
                step = it.step;
                step_row = it.step_row;
                return *this;
            }
        }

        mat_iterator<T> &operator+=(size_t move) noexcept(false) {
            for (size_t i = 0; i < move; i++) {
                ++*this;
            }
            return *this;
        }

        mat_iterator<T> &operator-=(size_t move) noexcept(false) {
            for (size_t i = 0; i < move; i++) {
                --*this;
            }
            return *this;
        }

        long long operator-(const mat_iterator<T> &m1) {
            long long row_dif = cur_row - m1.cur_row;
            long long col_dif = cur_col - m1.cur_col;
            long long chan_dif = cur_chan - m1.cur_chan;

            return (chan_dif * row_range * col_range) + (row_dif * col_range) + col_dif;
        }

        /*
         * 判断运算符
         */
        bool operator==(const mat_iterator<T> &it) const noexcept {
            return current == it.current;
        }

        bool operator!=(const mat_iterator<T> &it) const noexcept {
            return current != it.current;
        }

        /*
         * 解引用运算符
         */
        T &operator*() {
            return *current;
        }

    };


    typedef matrix<unsigned char> matUC8;
    typedef matrix<short> matI16;
    typedef matrix<int> matI32;
    typedef matrix<float> matF32;
    typedef matrix<double> matF64;

    /*
         * 1、构造器与析构函数
         */
    //不允许隐式转换
    template<typename T>
    matrix<T>::matrix(size_t rows, size_t cols, size_t channels) noexcept(true):
            rows(rows), cols(cols),
            channels(channels), step(this->cols),
            step_row(this->rows), shift(0),
            channel_range(0, (long long) channels) {
        for (int i = 0; i < channels; i++) {
            std::shared_ptr<T> tmp(new(std::nothrow) T[rows * cols], self_def_deleter);
            if (tmp.get() == nullptr) {
                cerr << "Warning: failed to allocate memory, thus the " << i
                     << "th channel pointer of matrix is nullptr!" << endl;

            }
            this->data.push_back(tmp);
        }
        mat_iterator<T> it1(*this, this->data[channel_range.getLB()].get() + shift, 0, 0, channel_range.getLB());
        mat_iterator<T> it2(*this, this->data[channel_range.getUB() - 1].get() + shift + (rows - 1) * step + cols,
                            rows - 1, cols - 1, channel_range.getUB() - 1);
        begin = it1;
        cur = it1;
        end = it2;
    }

    //使用已有数据进行初始化
    template<typename T>
    matrix<T>::matrix(size_t rows, size_t cols, size_t channels, const T *data) noexcept(true):
            rows(rows), cols(cols),
            channels(channels), step(this->cols),
            step_row(this->rows), shift(0),
            channel_range(0, (long long) channels) {
        if (data == nullptr) {
            cerr << "Warning: data == nullptr, this matrix's data was allocated "
                 << rows << " * " << cols << " * " << channels << " * "
                 << sizeof(T) << " bytes." << endl;
            for (size_t i = 0; i < channels; i++) {
                std::shared_ptr<T> tmp(new(std::nothrow) T[rows * cols], self_def_deleter);
                if (tmp.get() == nullptr) {
                    cerr << "Warning: failed to allocate memory, thus the " << i <<
                         "th channel pointer of matrix is nullptr!" << endl;
                }
                this->data.push_back(tmp);
            }
        } else {
            for (size_t j = 0; j < channels; j++) {
                std::shared_ptr<T> tmp(new(std::nothrow) T[rows * cols], self_def_deleter);
                if (tmp.get() == nullptr) {
                    cerr << "Warning: failed to allocate memory, thus the " << j <<
                         "th channel pointer of matrix is nullptr!" << endl;
                }
                this->data.push_back(tmp);
                for (size_t i = 0; i < rows * cols; i++) {
                    this->data[j].get()[i] = data[j * rows * cols + i];
                }
            }
        }
        mat_iterator<T> it1(*this, this->data[channel_range.getLB()].get() + shift, 0, 0, channel_range.getLB());
        mat_iterator<T> it2(*this, this->data[channel_range.getUB() - 1].get() + shift + (rows - 1) * step + cols,
                            rows - 1, cols - 1, channel_range.getUB() - 1);
        begin = it1;
        cur = it1;
        end = it2;
    }

    //使用初始化列表进行初始化
    template<typename T>
    matrix<T>::matrix(size_t rows, size_t cols, size_t channels, const std::initializer_list<T> &list) noexcept(true):
            rows(rows), cols(cols),
            channels(channels), step(this->cols),
            step_row(this->rows), shift(0),
            channel_range(0, (long long) channels) {
        for (size_t j = 0; j < channels; j++) {
            std::shared_ptr<T> tmp(new(std::nothrow) T[rows * cols], self_def_deleter);
            if (tmp.get() == nullptr) {
                cerr << "Warning: failed to allocate memory, thus the " << j <<
                     "th channel pointer of matrix is nullptr!" << endl;
            }
            this->data.push_back(tmp);
        }
        int i = 0, j = 0;
        int k = 0;
        for (auto it = list.begin(); it != list.end() && k != rows * cols * channels; it++, k++) {
            data[i].get()[j] = *it;
            if (i == channels - 1) {
                i = 0;
                j++;
            } else i++;
        }

        mat_iterator<T> it1(*this, data[channel_range.getLB()].get() + shift, 0, 0, channel_range.getLB());
        mat_iterator<T> it2(*this, data[channel_range.getUB() - 1].get() + shift + (rows - 1) * step + cols,
                            rows - 1, cols - 1, channel_range.getUB() - 1);
        begin = it1;
        cur = it1;
        end = it2;
    }

    //全部设为一个值
    template<typename T>
    matrix<T>::matrix(size_t rows, size_t cols, size_t channels, T t) noexcept(true):
            rows(rows), cols(cols),
            channels(channels), step(this->cols),
            step_row(this->rows), shift(0),
            channel_range(0, (long long) channels) {
        for (size_t j = 0; j < channels; j++) {
            std::shared_ptr<T> tmp(new(std::nothrow) T[rows * cols], self_def_deleter);
            if (tmp.get() == nullptr) {
                cerr << "Warning: failed to allocate memory, thus the " << j <<
                     "th channel pointer of matrix is nullptr!" << endl;
            }
            this->data.push_back(tmp);
            for (int i = 0; i < rows * cols; i++) {
                data[j].get()[i] = t;
            }
        }
        mat_iterator<T> it1(*this, data[channel_range.getLB()].get() + shift, 0, 0, channel_range.getLB());
        T *pt = data[channel_range.getUB() - 1].get() + shift + (rows - 1) * step + cols;
        mat_iterator<T> it2(*this, pt, rows - 1, cols - 1, channel_range.getUB() - 1);
        begin = it1;
        cur = it1;
        end = it2;
    }

    //创建特殊矩阵
    template<typename T>
    matrix<T>::matrix(size_t rows, size_t channels, mat_type type) noexcept(true):
            rows(rows), cols(rows),
            channels(channels), step(this->cols),
            step_row(this->rows), shift(0),
            channel_range(0, (long long) channels) {
        switch (type) {
            case identity_mat: {
                for (size_t i = 0; i < channels; i++) {
                    std::shared_ptr<T> tmp(new(std::nothrow) T[rows * rows], self_def_deleter);
                    if (tmp.get() == nullptr) {
                        cerr << "Warning: failed to allocate memory, thus the " << i <<
                             "th channel pointer of matrix is nullptr!" << endl;
                    }
                    this->data.push_back(tmp);
                    T *pt = data[i].get();
                    memset(pt, 0, rows * rows * sizeof(T));
                    for (int j = 0; j < rows; j++) {
                        pt[j * cols + j] = 1;
                    }
                }
                break;
            }
            default: {
                cerr << "Warning: Unsupported mode to create such matrix! This matrix's data was allocated "
                     << 1 << " * " << 1 << " * " << channels << " * "
                     << sizeof(T) << " bytes."
                     << "remark that the elements in the data is not initialized" << endl;
                for (size_t i = 0; i < channels; i++) {
                    std::shared_ptr<T> tmp(new(std::nothrow) T[1 * 1], self_def_deleter);
                    if (tmp.get() == nullptr) {
                        cerr << "Warning: failed to allocate memory, thus the " << i <<
                             "th channel pointer of matrix is nullptr!" << endl;
                    }
                    this->data.push_back(tmp);
                }
            }
        }
        mat_iterator<T> it1(*this, data[channel_range.getLB()].get() + shift, 0, 0, channel_range.getLB());
        mat_iterator<T> it2(*this, data[channel_range.getUB() - 1].get() + shift + (rows - 1) * step + cols,
                            rows - 1, cols - 1, channel_range.getUB() - 1);
        begin = it1;
        cur = it1;
        end = it2;
    }

    //拷贝构造器
    template<typename T>
    matrix<T>::matrix(const matrix &m) noexcept(true):
            rows(m.rows), cols(m.cols),
            channels(m.channels), step(m.step),
            step_row(m.step_row), shift(m.shift),
            channel_range(m.channel_range.getLB(), m.channel_range.getUB()),
            begin(m.begin), cur(m.cur), end(m.end) {
        data = m.data;
    }

    //ROI截取构造器
    template<typename T>
    matrix<T>::matrix(matrix &m, mat_range rowRange, mat_range colRange, mat_range chanRange) noexcept(false) {
#ifndef __NOEXCEPT_TRUE__
        if(rowRange.getLB() >= rowRange.getUB())
        {
            InvalidRangeException e(rowRange.getLB(), rowRange.getUB(), __WHERE__);
            throw e;
        }
        else if(colRange.getLB() >= colRange.getUB())
        {
            InvalidRangeException e(colRange.getLB(), colRange.getUB(), __WHERE__);
            throw e;
        }
        else if(chanRange.getLB() >= chanRange.getUB())
        {
            InvalidRangeException e(chanRange.getLB(), chanRange.getUB(), __WHERE__);
            throw e;
        }
        else if(rowRange.getLB() < 0)
        {
            IndexOutOfBoundException e(rowRange.getLB(), 0, __WHERE__);
            throw e;
        }
        else if(rowRange.getUB() > rows)
        {
            IndexOutOfBoundException e(rowRange.getLB(), rows, __WHERE__);
            throw e;
        }
        else if(colRange.getLB() < 0)
        {
            IndexOutOfBoundException e(colRange.getLB(), 0, __WHERE__);
            throw e;
        }
        else if(colRange.getUB() > cols)
        {
            IndexOutOfBoundException e(colRange.getLB(), cols, __WHERE__);
            throw e;
        }
        else if(channel_range.getLB() < m.channel_range.getLB())
        {
            IndexOutOfBoundException e(channel_range.getLB(), m.channel_range.getLB(), __WHERE__);
            throw e;
        }
        else if(channel_range.getUB() > m.channel_range.getUB())
        {
            IndexOutOfBoundException e(channel_range.getUB(), m.channel_range.getUB(), __WHERE__);
            throw e;
        }
#else
        if (rowRange.getLB() >= rowRange.getUB()) {
            InvalidRangeException e(rowRange.getLB(), rowRange.getUB(), __WHERE__);
            cout << endl;
            cerr << e.what() << endl;
            rowRange = mat_range(0, 0);
        } else if (colRange.getLB() >= colRange.getUB()) {
            InvalidRangeException e(colRange.getLB(), colRange.getUB(), __WHERE__);
            cout << endl;
            cerr << e.what() << endl;
            colRange = mat_range(0, 0);
        } else if (chanRange.getLB() >= chanRange.getUB()) {
            InvalidRangeException e(chanRange.getLB(), chanRange.getUB(), __WHERE__);
            cout << endl;
            cerr << e.what() << endl;
            chanRange = mat_range(0, 0);
        } else if (rowRange.getLB() < 0) {
            IndexOutOfBoundException e(rowRange.getLB(), 0, __WHERE__);
            cout << endl;
            cerr << e.what() << endl;
            rowRange = mat_range(0, 0);
        } else if (rowRange.getUB() > m.rows) {
            IndexOutOfBoundException e(rowRange.getUB(), rows, __WHERE__);
            cout << endl;
            cerr << e.what() << endl;
            rowRange = mat_range(0, 0);
        } else if (colRange.getLB() < 0) {
            IndexOutOfBoundException e(colRange.getLB(), 0, __WHERE__);
            cout << endl;
            cerr << e.what() << endl;
            colRange = mat_range(0, 0);
        } else if (colRange.getUB() > m.cols) {
            IndexOutOfBoundException e(colRange.getLB(), cols, __WHERE__);
            cout << endl;
            cerr << e.what() << endl;
            colRange = mat_range(0, 0);
        } else if (chanRange.getLB() < m.channel_range.getLB()) {
            IndexOutOfBoundException e(chanRange.getLB(), m.channel_range.getLB(), __WHERE__);
            cout << endl;
            cerr << e.what() << endl;
            chanRange = mat_range(0, 0);
        } else if (chanRange.getUB() > m.channel_range.getUB()) {
            //cout<<m.channel_range.getUB()<<endl;
            IndexOutOfBoundException e(chanRange.getUB(), m.channel_range.getUB(), __WHERE__);
            cout << endl;
            cerr << e.what() << endl;
            chanRange = mat_range(0, 0);
        }
#endif
        rows = rowRange.getUB() - rowRange.getLB();
        cols = colRange.getUB() - colRange.getLB();
        data = m.data;
        channel_range = chanRange;
        channels = chanRange.getUB() - chanRange.getLB();
        step = m.step;
        step_row = m.step_row;
        shift = m.shift + rowRange.getLB() * step + colRange.getLB();
        mat_iterator<T> it1(*this, data[channel_range.getLB()].get() + shift, 0, 0, channel_range.getLB());
        mat_iterator<T> it2(*this, data[channel_range.getLB()].get() + shift + (rows - 1) * step + cols,
                            0, 0, channel_range.getUB() - 1);
        begin = it1;
        cur = it1;
        end = it2;
    }


    //增加一层矩阵
    template<typename T>
    void matrix<T>::add_channel(std::initializer_list<T> list) noexcept {
        std::shared_ptr<T> tmp(new(std::nothrow) T[rows * cols], self_def_deleter);
        if (tmp.get() == nullptr) {
            cerr << "Warning: failed to allocate memory, thus the " << channels <<
                 "th channel pointer of matrix is nullptr!" << endl;
        }
        this->data.push_back(tmp);
        channels++;
        channel_range.ub++;
        int k = 0;
        for (auto it = list.begin(); it != list.end() && k != rows * cols; it++, k++) {
            data[channels - 1].get()[k] = *it;
        }
    }

    //删除指定层的矩阵
    template<typename T>
    void matrix<T>::erase_channel(size_t pos) noexcept {
        if (pos < channels) {
            data.erase(data.begin() + channel_range.getLB() + pos);
            channels--;
            channel_range.ub--;
        }
    }

    //ROI矩阵变大变小
    template<typename T>
    void matrix<T>::reshape(long long delta_row, long long delta_col, long long delta_chan) noexcept(false) {
#ifdef __NOEXCEPT_TRUE__
        //最少size为0
        if (cur.row_range + 2 * delta_row < 0)
            delta_row = cur.row_range / 2;
        if (cur.row_range.getUB() + delta_row > step_row) {
            //end指针移动
            end.current += (step_row - cur.row_range.getUB()) * step;
            end.cur_row = step_row;
            //范围变化
            cur.row_range.ub = step_row;
            begin.row_range.ub = step_row;
            end.row_range.ub = step_row;
        } else {
            //end指针移动
            end.current += delta_row * step;
            end.cur_row += delta_row;
            //范围变化
            cur.row_range.ub += delta_row;
            begin.row_range.ub += delta_row;
            end.row_range.ub += delta_row;
        }
        if (cur.row_range.getLB() - (long long) delta_row < 0) {
            //begin指针移动
            begin.current -= cur.row_range.getLB() * step;
            begin.cur_row -= cur.row_range.getLB();
            //迭代器的shift改变
            cur.shift -= cur.row_range.getLB() * step;
            begin.shift -= begin.row_range.getLB() * step;
            end.shift -= end.row_range.getLB() * step;
            //矩阵的shift改变
            shift -= end.row_range.getLB() * step;
            //范围变化
            cur.row_range.lb = 0;
            begin.row_range.lb = 0;
            end.row_range.lb = 0;
        } else {
            //begin指针移动
            begin.current -= delta_row * step;
            begin.cur_row -= delta_row;
            //迭代器的shift改变
            cur.shift -= delta_row * step;
            begin.shift -= delta_row * step;
            end.shift -= delta_row * step;
            //矩阵的shift改变
            shift -= delta_row * step;
            //范围变化
            begin.row_range.lb -= delta_row;
            end.row_range.lb -= delta_row;
            cur.row_range.lb -= delta_row;
        }
        rows = begin.row_range;

        //最少size为0
        if (cur.col_range + 2 * delta_col < 0)
            delta_col = cur.col_range / 2;
        if (cur.col_range.getUB() + delta_col > step) {
            //end指针移动
            end.current += step - cur.col_range.getUB();
            //范围变化
            cur.col_range.ub = step;
            begin.col_range.ub = step;
            end.col_range.ub = step;
        } else {
            //end指针移动
            end.current += delta_col;
            //范围变化
            cur.col_range.ub += delta_col;
            begin.col_range.ub += delta_col;
            end.col_range.ub += delta_col;
        }
        if (cur.col_range.getLB() - (long long) delta_col < 0) {
            //begin指针移动
            begin.current -= cur.col_range.getLB();
            //迭代器的shift改变
            cur.shift -= cur.col_range.getLB();
            begin.shift -= begin.col_range.getLB();
            end.shift -= end.col_range.getLB();
            //矩阵的shift改变
            shift -= end.col_range.getLB();
            //范围变化
            cur.col_range.lb = 0;
            begin.col_range.lb = 0;
            end.col_range.lb = 0;
        } else {
            //begin指针移动
            begin.current -= delta_col;
            //迭代器的shift改变
            cur.shift -= delta_col;
            begin.shift -= delta_col;
            end.shift -= delta_col;
            //矩阵的shift改变
            shift -= delta_col;
            //范围变化
            begin.col_range.lb -= delta_col;
            end.col_range.lb -= delta_col;
            cur.col_range.lb -= delta_col;
        }
        cols = begin.col_range;

        //最少size为0
        if (cur.channel_range + 2 * delta_chan < 0)
            delta_chan = cur.channel_range / 2;
        if (cur.channel_range.getUB() + delta_chan > cur.begin.size()) {
            //end指针移动
            end.current = cur.begin[cur.begin.size() - 1].get()
                          + cur.shift + (cur.row_range - 1) * step + cur.col_range;
            //范围变化
            cur.channel_range.ub = cur.begin.size();
            begin.channel_range.ub = cur.begin.size();
            end.channel_range.ub = cur.begin.size();
            channel_range.ub = cur.begin.size();
        } else {
            //end指针移动
            end.current = cur.begin[cur.channel_range.getUB() + delta_chan - 1].get()
                          + cur.shift + (cur.row_range - 1) * step + cur.col_range;
            //范围变化
            cur.channel_range.ub += delta_chan;
            begin.channel_range.ub += delta_chan;
            end.channel_range.ub += delta_chan;
            channel_range.ub += delta_chan;
        }
        if (cur.channel_range.getLB() - (long long) delta_chan < 0) {
            //begin指针移动
            begin.current = cur.begin[0].get() + cur.shift;
            //范围变化
            cur.channel_range.lb = 0;
            begin.channel_range.lb = 0;
            end.channel_range.lb = 0;
            channel_range.lb = 0;
        } else {
            //begin指针移动
            begin.current = cur.begin[cur.channel_range.lb - delta_chan].get() + cur.shift;
            //范围变化
            begin.channel_range.lb -= delta_chan;
            end.channel_range.lb -= delta_chan;
            cur.channel_range.lb -= delta_chan;
            channel_range.lb -= delta_chan;
        }
        channels = begin.channel_range;
#else
        //最少size为0
        if (cur.row_range + 2 * delta_row < 0) {
            InvalidRangeException e(cur.row_range.getUB() + delta_row, cur.row_range.getLB() - delta_row, __WHERE__);
            throw e;
        }
        if (cur.row_range.getUB() + delta_row > step_row) {
            IndexOutOfBoundException e(cur.row_range.getUB() + delta_row, step_row, __WHERE__);
            throw e;
        }
        if (cur.row_range.getLB() - (long long) delta_row < 0) {
            IndexOutOfBoundException e(cur.row_range.getLB() - (long long) delta_row, 0, __WHERE__);
            throw e;
        }
        //最少size为0
        if (cur.col_range + 2 * delta_col < 0) {
            InvalidRangeException e(cur.col_range.getUB() + delta_col, cur.col_range.getLB() - delta_col, __WHERE__);
            throw e;
        }
        if (cur.col_range.getUB() + delta_col > step) {
            IndexOutOfBoundException e(cur.col_range.getUB() + delta_col, step, __WHERE__);
            throw e;
        }
        if (cur.col_range.getLB() - (long long) delta_col < 0) {
            IndexOutOfBoundException e(cur.col_range.getLB() - (long long) delta_col, 0, __WHERE__);
            throw e;
        }
        //最少size为0
        if (cur.channel_range + 2 * delta_chan < 0) {
            InvalidRangeException e(cur.channel_range.getUB() + delta_chan, cur.channel_range.getLB() - delta_chan,
                                    __WHERE__);
            throw e;
        }
        if (cur.channel_range.getUB() + delta_chan > cur.begin.size()) {
            IndexOutOfBoundException e(cur.channel_range.getLB() + delta_chan, cur.begin.size(), __WHERE__);
            throw e;
        }
        if (cur.channel_range.getLB() - (long long) delta_chan < 0) {
            IndexOutOfBoundException e(cur.channel_range.getLB() - (long long) delta_chan, 0, __WHERE__);
            throw e;
        }
        //end指针移动
        end.current += delta_row * step;
        //范围变化
        cur.row_range.ub += delta_row;
        begin.row_range.ub += delta_row;
        end.row_range.ub += delta_row;
        //begin指针移动
        begin.current -= delta_row * step;
        //迭代器的shift改变
        cur.shift -= delta_row * step;
        begin.shift -= delta_row * step;
        end.shift -= delta_row * step;
        //矩阵的shift改变
        shift -= delta_row * step;
        //范围变化
        begin.row_range.lb -= delta_row;
        end.row_range.lb -= delta_row;
        cur.row_range.lb -= delta_row;
        rows = begin.row_range;
        //end指针移动
        end.current += delta_col;
        //范围变化
        cur.col_range.ub += delta_col;
        begin.col_range.ub += delta_col;
        end.col_range.ub += delta_col;
        //begin指针移动
        begin.current -= delta_col;
        //迭代器的shift改变
        cur.shift -= delta_col;
        begin.shift -= delta_col;
        end.shift -= delta_col;
        //矩阵的shift改变
        shift -= delta_col;
        //范围变化
        begin.col_range.lb -= delta_col;
        end.col_range.lb -= delta_col;
        cur.col_range.lb -= delta_col;
        cols = begin.col_range;
        //end指针移动
        end.current = cur.begin[cur.channel_range.getUB() + delta_chan - 1].get()
                      + cur.shift + (cur.row_range - 1) * step + cur.col_range;
        //范围变化
        cur.channel_range.ub += delta_chan;
        begin.channel_range.ub += delta_chan;
        end.channel_range.ub += delta_chan;
        channel_range.ub += delta_chan;
        //begin指针移动
        begin.current = cur.begin[cur.channel_range.lb - delta_chan].get() + cur.shift;
        //范围变化
        begin.channel_range.lb -= delta_chan;
        end.channel_range.lb -= delta_chan;
        cur.channel_range.lb -= delta_chan;
        channel_range.lb -= delta_chan;
        channels = begin.channel_range;
#endif
    }

    //ROI矩阵上下左右前后的移动
    template<typename T>
    void matrix<T>::move(long long delta_row, long long delta_col, long long delta_chan) noexcept(false) {
#ifdef __NOEXCEPT_TRUE__
        if (cur.row_range.getUB() + delta_row > step_row) {
            //迭代器改变
            end.current += (step_row - cur.row_range.getUB()) * step;
            end.cur_row = step_row;
            begin.current += (step_row - cur.row_range.getUB()) * step;
            begin.cur_row = step_row;
            cur.current += (step_row - cur.row_range.getUB()) * step;
            cur.cur_row = step_row;
            //row_range改变
            long long delta = step_row - cur.row_range.getUB();
            cur.row_range.ub = step_row;
            begin.row_range.ub = step_row;
            end.row_range.ub = step_row;
            cur.row_range.lb += delta;
            begin.row_range.lb += delta;
            end.row_range.lb += delta;
            //迭代器的shift改变
            cur.shift += delta * step;
            begin.shift += delta * step;
            end.shift += delta * step;
            //矩阵的shift改变
            shift += delta * step;
        } else if (cur.row_range.getLB() + delta_row < 0) {
            //迭代器改变
            begin.current -= cur.row_range.getLB() * step;
            begin.cur_row -= cur.row_range.getLB();
            end.current -= cur.row_range.getLB() * step;
            end.cur_row -= cur.row_range.getLB();
            cur.current -= cur.row_range.getLB() * step;
            cur.cur_row -= cur.row_range.getLB();
            //row_range改变
            long long delta = cur.row_range.getLB();
            cur.row_range.ub -= cur.row_range.getLB();
            begin.row_range.ub -= cur.row_range.getLB();
            end.row_range.ub -= cur.row_range.getLB();
            cur.row_range.lb = 0;
            begin.row_range.lb = 0;
            end.row_range.lb = 0;
            //迭代器的shift改变
            cur.shift -= delta * step;
            begin.shift -= delta * step;
            end.shift -= delta * step;
            //矩阵的shift改变
            shift -= delta * step;
        } else {
            //迭代器改变
            end.current += delta_row * step;
            end.cur_row += delta_row;
            begin.current += delta_row * step;
            begin.cur_row += delta_row;
            cur.current += delta_row * step;
            cur.cur_row += delta_row;
            //row_range改变
            cur.row_range.ub += delta_row;
            begin.row_range.ub += delta_row;
            end.row_range.ub += delta_row;
            cur.row_range.lb += delta_row;
            begin.row_range.lb += delta_row;
            end.row_range.lb += delta_row;
            //迭代器的shift改变
            cur.shift += delta_row * step;
            begin.shift += delta_row * step;
            end.shift += delta_row * step;
            //矩阵的shift改变
            shift += delta_row * step;
        }

        if (cur.col_range.getUB() + delta_col > step) {
            //迭代器改变
            end.current += (step - cur.col_range.getUB());
            end.cur_col = step;
            begin.current += (step - cur.col_range.getUB());
            begin.cur_col = step;
            cur.current += (step - cur.col_range.getUB());
            cur.cur_col = step;
            //row_range改变
            long long delta = step - cur.col_range.getUB();
            cur.col_range.ub = step;
            begin.col_range.ub = step;
            end.col_range.ub = step;
            cur.col_range.lb += delta;
            begin.col_range.lb += delta;
            end.col_range.lb += delta;
            //迭代器的shift改变
            cur.shift += delta;
            begin.shift += delta;
            end.shift += delta;
            //矩阵的shift改变
            shift += delta;
        } else if (cur.col_range.getLB() + delta_col < 0) {
            //迭代器改变
            begin.current -= cur.col_range.getLB();
            begin.cur_col -= cur.col_range.getLB();
            end.current -= cur.col_range.getLB();
            end.cur_col -= cur.col_range.getLB();
            cur.current -= cur.col_range.getLB();
            cur.cur_col -= cur.col_range.getLB();
            //col_range改变
            long long delta = cur.col_range.getLB();
            cur.col_range.ub -= cur.col_range.getLB();
            begin.col_range.ub -= cur.col_range.getLB();
            end.col_range.ub -= cur.col_range.getLB();
            cur.col_range.lb = 0;
            begin.col_range.lb = 0;
            end.col_range.lb = 0;
            //迭代器的shift改变
            cur.shift -= delta;
            begin.shift -= delta;
            end.shift -= delta;
            //矩阵的shift改变
            shift -= delta;
        } else {
            //迭代器改变
            end.current += delta_col;
            end.cur_col += delta_col;
            begin.current += delta_col;
            begin.cur_col += delta_col;
            cur.current += delta_col;
            cur.cur_col += delta_col;
            //col_range改变
            cur.col_range.ub += delta_col;
            begin.col_range.ub += delta_col;
            end.col_range.ub += delta_col;
            cur.col_range.lb += delta_col;
            begin.col_range.lb += delta_col;
            end.col_range.lb += delta_col;
            //迭代器的shift改变
            cur.shift += delta_col;
            begin.shift += delta_col;
            end.shift += delta_col;
            //矩阵的shift改变
            shift += delta_col;
        }

        long long begin_s = shift, cur_s = shift + (cur.current - begin.current), end_s =
                shift + (end.current - begin.current);
        if (cur.channel_range.getUB() + delta_chan > data.size()) {
            //迭代器改变
            long long delta = data.size() - cur.channel_range.getUB();
            end.current = data[end.cur_chan + delta].get() + end_s;
            begin.current = data[begin.cur_chan + delta].get() + begin_s;
            cur.current = data[cur.cur_chan + delta].get() + cur_s;
            //channel_range，cur_chan改变
            cur.channel_range.ub = data.size();
            begin.channel_range.ub = data.size();
            end.channel_range.ub = data.size();
            cur.channel_range.lb += delta;
            begin.channel_range.lb += delta;
            end.channel_range.lb += delta;
            cur.cur_chan += delta;
            begin.cur_chan += delta;
            end.cur_chan += delta;
            //矩阵的channel_range改变
            channel_range.ub = data.size();
            channel_range.lb += delta;
        } else if (cur.channel_range.getLB() + delta_chan < 0) {
            //迭代器改变
            long long delta = cur.channel_range.getLB();
            end.current = data[end.cur_chan + delta].get() + end_s;
            begin.current = data[begin.cur_chan + delta].get() + begin_s;
            cur.current = data[cur.cur_chan + delta].get() + cur_s;
            //channel_range改变
            cur.channel_range.ub -= delta;
            begin.channel_range.ub -= delta;
            end.channel_range.ub -= delta;
            cur.channel_range.lb = 0;
            begin.channel_range.lb = 0;
            end.channel_range.lb = 0;
            cur.cur_chan -= delta;
            begin.cur_chan -= delta;
            end.cur_chan -= delta;
            //矩阵的channel_range改变
            channel_range.ub -= delta;
            channel_range.lb = 0;
        } else {
            //迭代器改变
            end.current = data[end.cur_chan + delta_chan].get() + end_s;
            begin.current = data[begin.cur_chan + delta_chan].get() + begin_s;
            cur.current = data[cur.cur_chan + delta_chan].get() + cur_s;
            //channel_range，cur_chan改变
            cur.channel_range.ub += delta_chan;
            begin.channel_range.ub += delta_chan;
            end.channel_range.ub += delta_chan;
            cur.channel_range.lb += delta_chan;
            begin.channel_range.lb += delta_chan;
            end.channel_range.lb += delta_chan;
            cur.cur_chan += delta_chan;
            begin.cur_chan += delta_chan;
            end.cur_chan += delta_chan;
            //矩阵的channel_range改变
            channel_range.ub += delta_chan;
            channel_range.lb += delta_chan;
        }
#else
        if(cur.row_range.getUB() + delta_row > step_row) {
            IndexOutOfBoundException e(cur.row_range.getUB() + delta_row, step_row, __WHERE__);
            throw e;
        }
        if(cur.row_range.getLB() + delta_row < 0) {
            IndexOutOfBoundException e(cur.row_range.getLB() + delta_row, 0, __WHERE__);
            throw e;
        }
        if(cur.col_range.getUB() + delta_col > step) {
            IndexOutOfBoundException e(cur.col_range.getUB() + delta_col, step, __WHERE__);
            throw e;
        }
        else if(cur.col_range.getLB() + delta_col < 0) {
            IndexOutOfBoundException e(cur.col_range.getLB() + delta_col, 0, __WHERE__);
            throw e;
        }
        if(cur.channel_range.getUB() + delta_chan > data.size()) {
            IndexOutOfBoundException e(cur.channel_range.getUB() + delta_chan, data.size(), __WHERE__);
            throw e;
        }
        if(cur.channel_range.getLB() + delta_chan < 0) {
            IndexOutOfBoundException e(cur.channel_range.getLB() + delta_chan , 0, __WHERE__);
            throw e;
        }
        {
            //迭代器改变
            end.current += delta_row * step;
            end.cur_row += delta_row;
            begin.current += delta_row * step;
            begin.cur_row += delta_row;
            cur.current += delta_row * step;
            cur.cur_row += delta_row;
            //row_range改变
            cur.row_range.ub += delta_row;
            begin.row_range.ub += delta_row;
            end.row_range.ub += delta_row;
            cur.row_range.lb += delta_row;
            begin.row_range.lb += delta_row;
            end.row_range.lb += delta_row;
            //迭代器的shift改变
            cur.shift += delta_row * step;
            begin.shift += delta_row * step;
            end.shift += delta_row * step;
            //矩阵的shift改变
            shift += delta_row * step;

            //迭代器改变
            end.current += delta_col;
            end.cur_col += delta_col;
            begin.current += delta_col;
            begin.cur_col += delta_col;
            cur.current += delta_col;
            cur.cur_col += delta_col;
            //col_range改变
            cur.col_range.ub += delta_col;
            begin.col_range.ub += delta_col;
            end.col_range.ub += delta_col;
            cur.col_range.lb += delta_col;
            begin.col_range.lb += delta_col;
            end.col_range.lb += delta_col;
            //迭代器的shift改变
            cur.shift += delta_col;
            begin.shift += delta_col;
            end.shift += delta_col;
            //矩阵的shift改变
            shift += delta_col;

        long long begin_s = shift, cur_s = shift + (cur.current - begin.current), end_s = shift + (end.current - begin.current);

            //迭代器改变
            end.current = data[end.cur_chan + delta_chan].get() + end_s;
            begin.current = data[begin.cur_chan + delta_chan].get() + begin_s;
            cur.current = data[cur.cur_chan + delta_chan].get() + cur_s;
            //channel_range，cur_chan改变
            cur.channel_range.ub += delta_chan;
            begin.channel_range.ub += delta_chan;
            end.channel_range.ub += delta_chan;
            cur.channel_range.lb += delta_chan;
            begin.channel_range.lb += delta_chan;
            end.channel_range.lb += delta_chan;
            cur.cur_chan += delta_chan;
            begin.cur_chan += delta_chan;
            end.cur_chan += delta_chan;
            //矩阵的channel_range改变
            channel_range.ub += delta_chan;
            channel_range.lb += delta_chan;
        }
#endif
    }


    //判断是否是子矩阵
    template<typename T>
    bool matrix<T>::isSubMatrix(const matrix<T> &m1) const noexcept {
        return (m1.data[0].get() == data[0].get()) &&
               (m1.channel_range.getLB() <= channel_range.getLB()) &&
               (m1.channel_range.getUB() >= channel_range.getUB()) &&
               (m1.shift / m1.step <= shift / step) &&
               (m1.shift / m1.step + m1.rows >= shift / step + rows) &&
               (m1.shift % m1.step <= shift % step) &&
               (m1.shift % m1.step + m1.cols >= shift % step + cols);
    }

    //析构函数
    template<typename T>
    matrix<T>::~matrix() {
        for (size_t i = channel_range.getLB(); i < channel_range.getUB(); i++) {
            data[i].reset();
        }
    }


    //赋值操作
    template<typename T>
    matrix<T> &matrix<T>::operator=(const matrix<T> &matrix1) noexcept {
        if (this == &matrix1) {
            return *this;
        } else {
            this->rows = matrix1.rows;
            this->cols = matrix1.cols;
            this->data = matrix1.data;
            this->channels = matrix1.channels;
            this->step = matrix1.step;
            this->step_row = matrix1.step_row;
            this->shift = matrix1.shift;
            this->channel_range = matrix1.channel_range;
            this->begin = matrix1.begin;
            this->cur = matrix1.cur;
            this->end = matrix1.end;
            this->pm = matrix1.pm;
            return *this;
        }
    }

    //硬拷贝
    template<typename T>
    void matrix<T>::hardCopy(const matrix<T> &m1) noexcept {
        size_t final_row = rows < m1.rows ? rows : m1.rows,
                final_col = rows < m1.cols ? cols : m1.cols,
                final_chan = rows < m1.channels ? cols : m1.channels;
        for (size_t k = 0; k < final_chan; k++) {
            for (size_t i = 0; i < final_row; i++) {
                for (size_t j = 0; j < final_col; j++) {
                    (*this)(k, i, j) = m1(k, i, j);
                }
            }
        }
    }

    template<typename T>
    //访问第chan层的矩阵头指针，只能是右值
    T *matrix<T>::operator[](size_t chan) const noexcept(false) {
        if (chan >= channels) {
            IndexOutOfBoundException e(chan, this->rows, __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what() << endl;
            return data[channels - 1 + channel_range.getLB()].get() + shift;
#else
            throw e;
#endif
        } else {
            return data[chan + channel_range.getLB()].get() + shift;
        }
    }

    //访问第i层j行k列的元素,可赋值
    template<typename T>
    T &matrix<T>::operator()(size_t i, size_t j, size_t k) const noexcept(false) {
        if (i >= channels) {
            IndexOutOfBoundException e(i, channels, __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what() << endl;
            i = channels - 1;
#else
            throw e;
#endif
        } else if (j >= rows) {
            IndexOutOfBoundException e(j, this->rows, __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what() << endl;
            j = rows - 1;
#else
            throw e;
#endif
        } else if (k >= cols) {
            IndexOutOfBoundException e(k, this->cols, __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what() << endl;
            k = cols - 1;
#else
            throw e;
#endif
        }
        return (data[i + channel_range.getLB()].get())[shift + j * step + k];

    }

    //返回行数
    template<typename T>
    size_t matrix<T>::getRows() const noexcept {
        return rows;
    }

    //返回列数
    template<typename T>
    size_t matrix<T>::getCols() const noexcept {
        return cols;
    }

    //返回通道数
    template<typename T>
    size_t matrix<T>::getChannels() const noexcept {
        return channels;
    }

    //返回步长
    template<typename T>
    size_t matrix<T>::getStep() const noexcept {
        return step;
    }

    //返回父矩阵的行数
    template<typename T>
    size_t matrix<T>::getStepRows() const noexcept {
        return step_row;
    }

    //返回偏移量
    template<typename T>
    size_t matrix<T>::getShift() const noexcept {
        return shift;
    }


    //将第chan个通道全部设为一个值
    template<typename T>
    void matrix<T>::set_channel(size_t chan, T val) noexcept(false) {
        if (chan >= channels) {
            IndexOutOfBoundException e(chan, channels, __WHERE__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what() << endl;
            chan = channel_range.getUB() - 1;
#else
            throw e;
#endif
        }
        T *pt = data[chan].get() + shift;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                pt[i * step + j] = val;
            }
        }
    }

    //将第chan个通道第i行设为一个值
    template<typename T>
    void matrix<T>::set_row(size_t chan, size_t i, T val) noexcept(false) {
        if (chan >= channels) {
            IndexOutOfBoundException e(chan, channels, __WHERE__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what() << endl;
            chan = channel_range.getUB() - 1;
#else
            throw e;
#endif
        }
        if (i >= rows) {
            IndexOutOfBoundException e(i, rows, __WHERE__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what() << endl;
            i = rows - 1;
#else
            throw e;
#endif
        }
        T *pt = data[chan].get() + shift;
        for (int j = 0; j < cols; j++) {
            pt[i * step + j] = val;
        }
    }

    //将第chan个通道第j列设为一个值
    template<typename T>
    void matrix<T>::set_col(size_t chan, size_t j, T val) noexcept(false) {
        if (chan >= channels) {
            IndexOutOfBoundException e(chan, channels, __WHERE__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what() << endl;
            chan = channel_range.getUB() - 1;
#else
            throw e;
#endif
        }
        if (j >= cols) {
            IndexOutOfBoundException e(j, cols, __WHERE__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what() << endl;
            j = cols - 1;
#else
            throw e;
#endif
        }
        T *pt = data[chan].get() + shift;
        for (int i = 0; i < rows; i++) {
            pt[i * step + j] = val;
        }
    }

    /*
         * 加法
         */
    //矩阵和常数相加
    template<typename T>
    template<typename V>
    [[nodiscard]]  //开辟了新内存，不应舍弃
    matrix<T> matrix<T>::operator+(V val) const noexcept {
        matrix<T> m1(rows, cols, channel_range);
        for (size_t k = channel_range.getLB(); k < channel_range.getUB(); k++) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    m1.data[k - channel_range.getLB()].get()[m1.shift + i * m1.step + j] =
                            data[k].get()[shift + i * step + j] + static_cast<T>(val);
                }
            }
        }
        return m1;
    }

    //矩阵和矩阵相加
    template<typename T>
    template<class V>
    [[nodiscard]]  //开辟了新内存，不应舍弃
    matrix<T> matrix<T>::operator+(const matrix<V> &m) const noexcept(false) {
        size_t final_row = rows, final_col = cols, final_chan = channels;
        if (this->rows != m.rows) {
            OperationException e(
                    "Rows are unequal! rows1 = " + to_string(this->rows) + ", rows2 = " + to_string(m.rows),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_row = this->rows < m.rows ? this->rows : m.rows;
#else
            throw e;
#endif
        }
        if (this->cols != m.cols) {
            OperationException e(
                    "Cols are unequal! cols1 = " + to_string(this->cols) + ", cols2 = " + to_string(m.cols),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_col = this->cols < m.cols ? this->cols : m.cols;
#else
            throw e;
#endif
        }
        if (this->channels != m.channels) {
            OperationException e(
                    "Channels are unequal! channels1 = " + to_string(this->channels) + ", channels2 = " +
                    to_string(m.channels),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_chan = this->channels < m.channels ? this->channels : m.channels;
#else
            throw e;
#endif
        }
        {
            matrix<T> m1(final_row, final_col, final_chan);
            for (size_t k = 0; k < final_chan; k++) {
                for (size_t i = 0; i < final_row; i++) {
                    for (size_t j = 0; j < final_col; j++) {
                        m1.data[k].get()[m1.shift + i * m1.step + j] =
                                data[channel_range.getLB() + k].get()[shift + i * step + j]
                                + m.data[m.channel_range.getLB() + k].get()[m.shift + i * m.step + j];
                    }
                }
            }
            return m1;
        }
    }

    //运算符+=
    template<typename T>
    template<typename V>
    matrix<T> &matrix<T>::operator+=(V val) noexcept {
        for (size_t k = channel_range.getLB(); k < channel_range.getUB(); k++) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    data[k].get()[shift + i * step + j] += static_cast<T>(val);
                }
            }
        }
        return *this;
    }

    template<typename T>
    matrix<T> &matrix<T>::operator+=(const matrix<T> &m) noexcept(false) {
        size_t final_row = rows, final_col = cols, final_chan = channels;
        if (this->rows != m.rows) {
            OperationException e(
                    "Rows are unequal! rows1 = " + to_string(this->rows) + ", rows2 = " + to_string(m.rows),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_row = this->rows < m.rows ? this->rows : m.rows;
#else
            throw e;
#endif
        }
        if (this->cols != m.cols) {
            OperationException e(
                    "Cols are unequal! cols1 = " + to_string(this->cols) + ", cols2 = " + to_string(m.cols),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_col = this->cols < m.cols ? this->cols : m.cols;
#else
            throw e;
#endif
        }
        if (this->channels != m.channels) {
            OperationException e(
                    "Channels are unequal! channels1 = " + to_string(this->channels) + ", channels2 = " +
                    to_string(m.channels),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_chan = this->channels < m.channels ? this->channels : m.channels;
#else
            throw e;
#endif
        }
        for (size_t k = 0; k < final_chan; k++) {
            for (size_t i = 0; i < final_row; i++) {
                for (size_t j = 0; j < final_col; j++) {
                    data[k + channel_range.getLB()].get()[shift + i * step + j] +=
                            m.data[k + m.channel_range.getLB()].get()[m.shift + i * m.step + j];
                }
            }
        }
        return *this;
    }

    /*
     * 减法
     */
    //矩阵和常数相减
    template<typename T>
    template<typename V>
    [[nodiscard]]  //开辟了新内存，不应舍弃
    matrix<T> matrix<T>::operator-(V val) const noexcept {
        matrix<T> m1(rows, cols, channel_range);
        for (size_t k = channel_range.getLB(); k < channel_range.getUB(); k++) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    m1.data[k - channel_range.getLB()].get()[m1.shift + i * m1.step + j] =
                            data[k].get()[shift + i * step + j] - static_cast<T>(val);
                }
            }
        }
        return m1;
    }

    //矩阵和矩阵相减
    template<typename T>
    [[nodiscard]]  //开辟了新内存，不应舍弃
    matrix<T> matrix<T>::operator-(const matrix<T> &m) const noexcept(false) {
        size_t final_row = rows, final_col = cols, final_chan = channels;
        if (this->rows != m.rows) {
            OperationException e(
                    "Rows are unequal! rows1 = " + to_string(this->rows) + ", rows2 = " + to_string(m.rows),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_row = this->rows < m.rows ? this->rows : m.rows;
#else
            throw e;
#endif
        }
        if (this->cols != m.cols) {
            OperationException e(
                    "Cols are unequal! cols1 = " + to_string(this->cols) + ", cols2 = " + to_string(m.cols),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_col = this->cols < m.cols ? this->cols : m.cols;
#else
            throw e;
#endif
        }
        if (this->channels != m.channels) {
            OperationException e(
                    "Channels are unequal! channels1 = " + to_string(this->channels) + ", channels2 = " +
                    to_string(m.channels),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_chan = this->channels < m.channels ? this->channels : m.channels;
#else
            throw e;
#endif
        }
        {
            matrix<T> m1(final_row, final_col, final_chan);
            for (size_t k = 0; k < final_chan; k++) {
                for (size_t i = 0; i < final_row; i++) {
                    for (size_t j = 0; j < final_col; j++) {
                        m1.data[k].get()[m1.shift + i * m1.step + j] =
                                data[channel_range.getLB() + k].get()[shift + i * step + j]
                                - m.data[m.channel_range.getLB() + k].get()[m.shift + i * m.step + j];
                    }
                }
            }
            return m1;
        }
    }

    //运算符-=
    template<typename T>
    template<typename V>
    matrix<T> &matrix<T>::operator-=(V val) noexcept {
        for (size_t k = channel_range.getLB(); k < channel_range.getUB(); k++) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    data[k].get()[shift + i * step + j] -= static_cast<T>(val);
                }
            }
        }
        return *this;
    }

    template<typename T>
    matrix<T> &matrix<T>::operator-=(const matrix<T> &m) noexcept(false) {
        size_t final_row = rows, final_col = cols, final_chan = channels;
        if (this->rows != m.rows) {
            OperationException e(
                    "Rows are unequal! rows1 = " + to_string(this->rows) + ", rows2 = " + to_string(m.rows),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_row = this->rows < m.rows ? this->rows : m.rows;
#else
            throw e;
#endif
        }
        if (this->cols != m.cols) {
            OperationException e(
                    "Cols are unequal! cols1 = " + to_string(this->cols) + ", cols2 = " + to_string(m.cols),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_col = this->cols < m.cols ? this->cols : m.cols;
#else
            throw e;
#endif
        }
        if (this->channels != m.channels) {
            OperationException e(
                    "Channels are unequal! channels1 = " + to_string(this->channels) + ", channels2 = " +
                    to_string(m.channels),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_chan = this->channels < m.channels ? this->channels : m.channels;
#else
            throw e;
#endif
        }
        for (size_t k = 0; k < final_chan; k++) {
            for (size_t i = 0; i < final_row; i++) {
                for (size_t j = 0; j < final_col; j++) {
                    data[k + channel_range.getLB()].get()[shift + i * step + j] -=
                            m.data[k + m.channel_range.getLB()].get()[m.shift + i * m.step + j];
                }
            }
        }
        return *this;
    }

    /*
         * 乘法
         */
    //矩阵和常数相乘
    template<typename T>
    template<typename V>
    [[nodiscard]]  //开辟了新内存，不应舍弃
    matrix<T> matrix<T>::operator*(V val) const noexcept {
        matrix<T> m1(rows, cols, channel_range);
        for (size_t k = channel_range.getLB(); k < channel_range.getUB(); k++) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    m1.data[k - channel_range.getLB()].get()[m1.shift + i * m1.step + j] =
                            data[k].get()[shift + i * step + j] * static_cast<T>(val);
                }
            }
        }
        return m1;
    }

    //矩阵和矩阵相乘
    template<typename T>
    [[nodiscard]]  //开辟了新内存，不应舍弃
    matrix<T> matrix<T>::operator*(const matrix<T> &m) const noexcept(false) {
        size_t col_row = cols, final_chan = channels;
        if (this->channels != m.channels) {
            OperationException e(
                    "Channels are unequal! channels1 = " + to_string(this->channels) + ", channels2 = " +
                    to_string(m.channels),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_chan = this->channels < m.channels ? this->channels : m.channels;
#else
            throw e;
#endif
        }
        if (this->cols != m.rows) {
            OperationException e("Cols1 and rows2 are unequal! Cols1 = " + to_string(this->cols) + ", rows2 = " +
                                 to_string(m.rows),
                                 __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            col_row = this->cols < m.rows ? this->cols : m.rows;
#else
            throw e;
#endif
        } else {
            matrix<T> m1(rows, m.cols, final_chan);
            for (size_t p = 0; p < final_chan; p++) {
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < m.cols; j++) {
                        for (int k = 0; k < col_row; k++) {

                            m1.data[p].get()[m1.shift + i * m1.step + j] +=
                                    data[p + channel_range.getLB()].get()[shift + i * step + k] *
                                    m.data[p + m.channel_range.getLB()].get()[m.shift + k * m.step + j];
                        }
                    }
                }
            }
            return m1;
        }
    }

    /*
    * 按位运算符：在这里只重载 ^ & | ~ \<\< >> 这几个最基本的按位运算符的最基本形式
    */
    //矩阵与数之间的^ & |
    template<typename T>
    template<typename V>
    [[nodiscard]]  //开辟了新内存，不应舍弃
    matrix<T> matrix<T>::operator^(V val) const noexcept {
        matrix<T> m1(rows, cols, channels);
        for (size_t k = 0; k < channels; k++) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    m1.data[channel_range.getLB()].get()[m1.shift + i * m1.step + j] =
                            data[k + channel_range.getLB()].get()[shift + i * step + j] ^ static_cast<T>(val);
                }
            }
        }
        return m1;
    }

    template<typename T>
    template<typename V>
    [[nodiscard]]  //开辟了新内存，不应舍弃
    matrix<T> matrix<T>::operator&(V val) const noexcept {
        matrix<T> m1(rows, cols, channels);
        for (size_t k = 0; k < channels; k++) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; ++j) {
                    m1.data[channel_range.getLB()].get()[m1.shift + i * m1.step + j] =
                            data[k + channel_range.getLB()].get()[shift + i * step + j] & static_cast<T>(val);
                }
            }
        }
        return m1;
    }

    template<typename T>
    template<typename V>
    [[nodiscard]]  //开辟了新内存，不应舍弃
    matrix<T> matrix<T>::operator|(V val) const noexcept {
        matrix<T> m1(rows, cols, channels);
        for (size_t k = 0; k < channels; k++) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; ++j) {
                    m1.data[channel_range.getLB()].get()[m1.step + i * m1.step + j] =
                            data[k + channel_range.getLB()].get()[step + i * step + j] | static_cast<T>(val);
                }
            }
        }
        return m1;
    }

    //矩阵与矩阵之间的按位^ & |
    template<typename T>
    [[nodiscard]]  //开辟了新内存，不应舍弃
    matrix<T> matrix<T>::operator^(const matrix<T> &m) const noexcept(false) {
        size_t final_row = rows, final_col = cols, final_chan = channels;
        if (this->rows != m.rows) {
            OperationException e(
                    "Rows are unequal! rows1 = " + to_string(this->rows) + ", rows2 = " + to_string(m.rows),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_row = this->rows < m.rows ? this->rows : m.rows;
#else
            throw e;
#endif
        }
        if (this->cols != m.cols) {
            OperationException e(
                    "Cols are unequal! cols1 = " + to_string(this->cols) + ", cols2 = " + to_string(m.cols),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_col = this->cols < m.cols ? this->cols : m.cols;
#else
            throw e;
#endif
        }
        if (this->channels != m.channels) {
            OperationException e(
                    "Channels are unequal! channels1 = " + to_string(this->channels) + ", channels2 = " +
                    to_string(m.channels),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_chan = this->channels < m.channels ? this->channels : m.channels;
#else
            throw e;
#endif
        }
        matrix<T> m1(final_row, final_col, final_chan);
        for (size_t k = 0; k < channels; k++) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    m1.data[channel_range.getLB()].get()[m1.shift + i * m1.step + j] =
                            data[k + channel_range.getLB()].get()[shift + i * step + j] ^
                            m.data[k + m.channel_range.getLB()].get()[m.shift + i * m.step + j];
                }
            }
        }
        return m1;
    }

    template<typename T>
    [[nodiscard]]  //开辟了新内存，不应舍弃
    matrix<T> matrix<T>::operator&(const matrix<T> &m) const noexcept(false) {
        size_t final_row = rows, final_col = cols, final_chan = channels;
        if (this->rows != m.rows) {
            OperationException e(
                    "Rows are unequal! rows1 = " + to_string(this->rows) + ", rows2 = " + to_string(m.rows),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_row = this->rows < m.rows ? this->rows : m.rows;
#else
            throw e;
#endif
        }
        if (this->cols != m.cols) {
            OperationException e(
                    "Cols are unequal! cols1 = " + to_string(this->cols) + ", cols2 = " + to_string(m.cols),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_col = this->cols < m.cols ? this->cols : m.cols;
#else
            throw e;
#endif
        }
        if (this->channels != m.channels) {
            OperationException e(
                    "Channels are unequal! channels1 = " + to_string(this->channels) + ", channels2 = " +
                    to_string(m.channels),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_chan = this->channels < m.channels ? this->channels : m.channels;
#else
            throw e;
#endif
        }
        matrix<T> m1(final_row, final_col, final_chan);
        for (size_t k = 0; k < channels; k++) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; ++j) {
                    m1.data[channel_range.getLB()].get()[m1.shift + i * m1.step + j] =
                            data[k + channel_range.getLB()].get()[shift + i * step + j] &
                            m.data[k + m.channel_range.getLB()].get()[m.shift + i * m.step + j];
                }
            }
        }
        return m1;

    }

    template<typename T>
    [[nodiscard]]  //开辟了新内存，不应舍弃
    matrix<T> matrix<T>::operator|(const matrix<T> &m) const noexcept(false) {
        size_t final_row = rows, final_col = cols, final_chan = channels;
        if (this->rows != m.rows) {
            OperationException e(
                    "Rows are unequal! rows1 = " + to_string(this->rows) + ", rows2 = " + to_string(m.rows),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_row = this->rows < m.rows ? this->rows : m.rows;
#else
            throw e;
#endif
        }
        if (this->cols != m.cols) {
            OperationException e(
                    "Cols are unequal! cols1 = " + to_string(this->cols) + ", cols2 = " + to_string(m.cols),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_col = this->cols < m.cols ? this->cols : m.cols;
#else
            throw e;
#endif
        }
        if (this->channels != m.channels) {
            OperationException e(
                    "Channels are unequal! channels1 = " + to_string(this->channels) + ", channels2 = " +
                    to_string(m.channels),
                    __FILE__, __LINE__, __func__);
#ifdef __NOEXCEPT_TRUE__
            cout << endl;
            cerr << e.what();
            final_chan = this->channels < m.channels ? this->channels : m.channels;
#else
            throw e;
#endif
        }
        matrix<T> m1(final_row, final_col, final_chan);
        for (size_t k = 0; k < channels; k++) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; ++j) {
                    m1.data[channel_range.getLB()].get()[m1.shift + i * m1.step + j] =
                            data[k + channel_range.getLB()].get()[shift + i * step + j] |
                            m.data[k + m.channel_range.getLB()].get()[m.shift + i * m.step + j];
                }
            }
        }
        return m1;
    }

    /*
         * 比较运算符的重载
         */
    //==: 矩阵规模和每个矩阵元素必须相等
    template<typename T>
    bool matrix<T>::operator==(const matrix<T> &m) const noexcept(false) {
        if (this->rows != m.rows) {
            return false;
        }
        if (this->cols != m.cols) {
            return false;
        }
        if (this->channels != m.channels) {
            return false;
        }
        {
            bool result = true;
            for (size_t k = 0; k < channels; k++) {
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < cols; j++) {
                        result &= (bool) (data[k + channel_range.getLB()].get()[shift + i * step + j] ==
                                          m.data[k + m.channel_range.getLB()].get()[m.shift + i * m.step + j]);
                    }
                }
            }
            return result;
        }
    }

    //\!=: 规模不同或至少有一个元素不相等
    template<typename T>
    bool matrix<T>::operator!=(const matrix<T> &m) const noexcept(false) {
        return !(*this == m);
    }

    // >=: 规模相同的情况下，每个元素都要大于等于， >: 每个元素都要大于
    template<typename T>
    bool matrix<T>::operator>=(const matrix<T> &m) const noexcept(false) {
        if (this->rows != m.rows) {
            return false;
        }
        if (this->cols != m.cols) {
            return false;
        }
        {
            bool result = true;
            for (size_t k = 0; k < channels; k++) {
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < cols; j++) {
                        result &= (bool) (data[k + channel_range.getLB()].get()[shift + i * step + j] >=
                                          m.data[k + m.channel_range.getLB()].get()[m.shift + i * m.step + j]);
                    }
                }
            }
            return result;
        }
    }

    template<typename T>
    bool matrix<T>::operator>(const matrix<T> &m) const noexcept(false) {
        if (this->rows != m.rows) {
            return false;
        }
        if (this->cols != m.cols) {
            return false;
        }
        {
            bool result = true;
            for (size_t k = 0; k < channels; k++) {
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < cols; j++) {
                        result &= (bool) (data[k + channel_range.getLB()].get()[shift + i * step + j] >
                                          m.data[k + m.channel_range.getLB()].get()[m.shift + i * m.step + j]);
                    }
                }
            }
            return result;
        }
    }

    // <=: 每个元素都要大于等于， <: 每个元素都要大于
    template<typename T>
    bool matrix<T>::operator<=(const matrix<T> &m) const noexcept(false) {
        if (this->rows != m.rows) {
            return false;
        }
        if (this->cols != m.cols) {
            return false;
        }
        {
            bool result = true;
            for (size_t k = 0; k < channels; k++) {
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < cols; j++) {
                        result &= (bool) (data[k + channel_range.getLB()].get()[shift + i * step + j] <=
                                          m.data[k + m.channel_range.getLB()].get()[m.shift + i * m.step + j]);
                    }
                }
            }
            return result;
        }
    }

    template<typename T>
    bool matrix<T>::operator<(const matrix<T> &m) const noexcept(false) {
        if (this->rows != m.rows) {
            return false;
        }
        if (this->cols != m.cols) {
            return false;
        }
        {
            bool result = true;
            for (size_t k = 0; k < channels; k++) {
                for (int i = 0; i < rows; i++) {
                    for (int j = 0; j < cols; j++) {
                        result &= (bool) (data[k + channel_range.getLB()].get()[shift + i * step + j] <
                                          m.data[k + m.channel_range.getLB()].get()[m.shift + i * m.step + j]);
                    }
                }
            }
            return result;

        }
    }

    //矩阵自身移位操作
    //左移
    template<typename T>
    template<typename V>
    matrix<T> & matrix<T>::operator<<=(V val) {
        for(size_t k = 0; k<channels; k++){
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    data[k + channel_range.getLB()].get()[shift + i * step + j] <<= val;
                }
            }
        }
        return *this;
    }

    template<typename T>
    template<typename V>
    matrix<T> matrix<T>::operator<<(V val) {
        matrix<T> m1(rows, cols, channels);
        m1.hardCopy(*this);
        for(size_t k = 0; k<channels; k++){
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    m1.data[k].get()[m1.shift + i * m1.step + j] =
                            data[k + channel_range.getLB()].get()[shift + i * step + j] << val;
                }
            }
        }
        return m1;
    }

    //右移
    template<typename T>
    template<typename V>
    matrix<T> & matrix<T>::operator>>=(V val) {
        for(size_t k = 0; k<channels; k++){
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    data[k + channel_range.getLB()].get()[shift + i * step + j] >>= val;
                }
            }
        }
        return *this;
    }

    template<typename T>
    template<typename V>
    matrix<T> matrix<T>::operator>>(V val) {
        matrix<T> m1(rows, cols, channels);
        m1.hardCopy(*this);
        for(size_t k = 0; k<channels; k++){
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    m1.data[k].get()[m1.shift + i * m1.step + j] =
                            data[k + channel_range.getLB()].get()[shift + i * step + j] >> val;
                }
            }
        }
        return m1;
    }
}


#endif //CPPPROJECTS_MATRIX_HPP
