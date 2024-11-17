#include "MATexception.hpp"

namespace mat {
    /*
     * 1、规模异常
     */
    InvalidSizeException::InvalidSizeException(long long wrong_size, string file_name,
                                               int line_number, string function_name) :
            wrong_size(wrong_size), file_name(std::move(file_name)),
            line_number(line_number), function_name(std::move(function_name)) {
        message = "InvalidSizeException occurs in file: " + this->file_name +
                  ", in function: " + this->function_name +
                  ", in line " + to_string(this->line_number) + "\n" +
                  "The wrong size is " + to_string(this->wrong_size) + "\n";
    }

    [[nodiscard]]
    const char *InvalidSizeException::what() const noexcept {
        return message.c_str();
    }

    [[nodiscard]]
    long long InvalidSizeException::getWrongSize() const noexcept {
        return wrong_size;
    }

    [[nodiscard]]
    const string &InvalidSizeException::getFileName() const noexcept {
        return file_name;
    }

    [[nodiscard]]
    int InvalidSizeException::getLineNumber() const noexcept {
        return line_number;
    }

    [[nodiscard]]
    const string &InvalidSizeException::getFunctionName() const noexcept {
        return function_name;
    }

    /*
     * 2、空指针异常
     */

    NullPointerException::NullPointerException(string parameter_name, string file_name, int line_number,
                                               string function_name) :
            parameter_name(std::move(parameter_name)), file_name(std::move(file_name)),
            line_number(line_number), function_name(std::move(function_name)) {
        message = "NullPointerException occurs in file: " + this->file_name +
                  ", in function: " + this->function_name +
                  ", in line " + to_string(this->line_number) + "\n" +
                  "The parameter name is " + this->parameter_name + "\n";
    }

    //程序终止时会打印
    [[nodiscard]]
    const char *NullPointerException::what() const noexcept {
        return message.c_str();
    }

    //用户自己做异常处理的时候可以有选择性的打印自己关注的信息，而不用每次都打印what();
    [[nodiscard]]
    const string &NullPointerException::getParameterName() const noexcept {
        return parameter_name;
    }

    [[nodiscard]]
    const string &NullPointerException::getFileName() const noexcept {
        return file_name;
    }

    [[nodiscard]]
    int NullPointerException::getLineNumber() const noexcept {
        return line_number;
    }

    [[nodiscard]]
    const string &NullPointerException::getFunctionName() const noexcept {
        return function_name;
    }

    /*
     * 3、越界异常
     */

    //构造器
    IndexOutOfBoundException::IndexOutOfBoundException
            (long long index, long long bound, string file_name, int line_number, string function_name) :
            index(index), bound(bound), file_name(std::move(file_name)),
            line_number(line_number), function_name(std::move(function_name)) {
        message = "IndexOutOfBoundException occurs in file: " + this->file_name +
                  ", in function: " + this->function_name +
                  ", in line " + to_string(this->line_number) + "\n" +
                  "Index " + to_string(this->index) + " out of bound " + to_string(this->bound) + "\n";
    }

    //程序终止时会打印
    [[nodiscard]]
    const char *IndexOutOfBoundException::what() const noexcept {
        return message.c_str();
    }

    //用户自己做异常处理的时候可以有选择性的打印自己关注的信息，而不用每次都打印what();
    [[nodiscard]]
    long long int IndexOutOfBoundException::getIndex() const noexcept {
        return index;
    }

    [[nodiscard]]
    long long int IndexOutOfBoundException::getBound() const noexcept {
        return bound;
    }

    [[nodiscard]]
    const string &IndexOutOfBoundException::getFileName() const noexcept {
        return file_name;
    }

    [[nodiscard]]
    int IndexOutOfBoundException::getLineNumber() const noexcept {
        return line_number;
    }

    [[nodiscard]]
    const string &IndexOutOfBoundException::getFunctionName() const noexcept {
        return function_name;
    }

    /*
     * 4、运算异常：处理加减法中矩阵size不相等，乘法中第一个列数和第二个的行数不相等
     */
    OperationException::OperationException(string fatal_description,
                                           string file_name, int line_number, string function_name) :
            fatal_description(std::move(fatal_description)), file_name(std::move(file_name)),
            line_number(line_number), function_name(std::move(function_name)) {
        message = "IndexOutOfBoundException occurs in file: " + this->file_name +
                  ", in function: " + this->function_name +
                  ", in line " + to_string(this->line_number) + "\n" +
                  this->fatal_description + "\n";
    }

    //程序终止时会打印
    [[nodiscard]]
    const char *OperationException::what() const noexcept {
        return message.c_str();
    }

    //用户自己做异常处理的时候可以有选择性的打印自己关注的信息，而不用每次都打印what();
    [[nodiscard]]
    const string &OperationException::getFatalDescription() const noexcept {
        return fatal_description;
    }

    [[nodiscard]]
    const string &OperationException::getFileName() const noexcept {
        return file_name;
    }

    [[nodiscard]]
    int OperationException::getLineNumber() const noexcept {
        return line_number;
    }

    [[nodiscard]]
    const string &OperationException::getFunctionName() const noexcept {
        return function_name;
    }

    /*
     * 5、迭代器越界异常
     */
    IteratorOutOfBoundException::IteratorOutOfBoundException(string fatal_description,
                                                             string file_name, int line_number, string function_name) :
            fatal_description(std::move(fatal_description)), file_name(std::move(file_name)),
            line_number(line_number), function_name(std::move(function_name)) {
        message = "IndexOutOfBoundException occurs in file: " + this->file_name +
                  ", in function: " + this->function_name +
                  ", in line " + to_string(this->line_number) + "\n" +
                  this->fatal_description + "\n";
    }

    //程序终止时会打印
    [[nodiscard]]
    const char *IteratorOutOfBoundException::what() const noexcept {
        return message.c_str();
    }

    //用户自己做异常处理的时候可以有选择性的打印自己关注的信息，而不用每次都打印what();
    [[nodiscard]]
    const string &IteratorOutOfBoundException::getFatalDescription() const noexcept {
        return fatal_description;
    }

    [[nodiscard]]
    const string &IteratorOutOfBoundException::getFileName() const noexcept {
        return file_name;
    }

    [[nodiscard]]
    int IteratorOutOfBoundException::getLineNumber() const noexcept {
        return line_number;
    }

    [[nodiscard]]
    const string &IteratorOutOfBoundException::getFunctionName() const noexcept {
        return function_name;
    }

    /*
     * 6、非方阵导致的越界异常
     */

    NotASquareMatrixException::NotASquareMatrixException(size_t rows, size_t cols,
                                                         string file_name, int line_number, string function_name) :
            rows(rows), cols(cols), file_name(std::move(file_name)), line_number(line_number),
            function_name(std::move(function_name)) {
        message = "IndexOutOfBoundException occurs in file: " + this->file_name +
                  ", in function: " + this->function_name +
                  ", in line " + to_string(this->line_number) + "\n" +
                  "rows = " + to_string(rows) + ", cols = " + to_string(cols) + "\n";
    }

    //程序终止时会打印
    [[nodiscard]]
    const char *NotASquareMatrixException::what() const noexcept {
        return message.c_str();
    }

    //用户自己做异常处理的时候可以有选择性的打印自己关注的信息，而不用每次都打印what();
    [[nodiscard]]
    const size_t &NotASquareMatrixException::getRows() const noexcept {
        return rows;
    }

    [[nodiscard]]
    const size_t &NotASquareMatrixException::getCols() const noexcept {
        return cols;
    }

    [[nodiscard]]
    const string &NotASquareMatrixException::getFileName() const noexcept {
        return file_name;
    }

    [[nodiscard]]
    int NotASquareMatrixException::getLineNumber() const noexcept {
        return line_number;
    }

    [[nodiscard]]
    const string &NotASquareMatrixException::getFunctionName() const noexcept {
        return function_name;
    }

    /*
     * 7、范围不合法
     */

    InvalidRangeException::InvalidRangeException
    (size_t lb, size_t ub,
     string file_name, int line_number, string function_name):
            lb(lb), ub(ub), file_name(std::move(file_name)), line_number(line_number),
            function_name(std::move(function_name))
    {
        message = "IndexOutOfBoundException occurs in file: " + this->file_name +
                  ", in function: " + this->function_name +
                  ", in line " + to_string(this->line_number) + "\n" +
                  "lower bound = " + to_string(lb) + ", upper bound = " + to_string(ub) + "\n";
    }

        //程序终止时会打印
        [[nodiscard]]
        const char* InvalidRangeException::what () const noexcept
        {
            return message.c_str();
        }

        //用户自己做异常处理的时候可以有选择性的打印自己关注的信息，而不用每次都打印what();
        [[nodiscard]]
        const size_t& InvalidRangeException::getUb() const noexcept
        {
            return ub;
        }

        [[nodiscard]]
        const size_t& InvalidRangeException::getLb() const noexcept
        {
            return lb;
        }

        [[nodiscard]]
        const string&  InvalidRangeException::getFileName() const noexcept
        {
            return file_name;
        }

        [[nodiscard]]
        int InvalidRangeException::getLineNumber() const noexcept
        {
            return line_number;
        }

        [[nodiscard]]
        const string& InvalidRangeException::getFunctionName() const noexcept
        {
            return function_name;
        }


}