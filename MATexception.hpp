#ifndef CPPPROJECTS_MATEXCEPTION_HPP
#define CPPPROJECTS_MATEXCEPTION_HPP
#include <exception>
#include <string>
//#include <utility>
using std::exception;
using std::string;
using std::to_string;
namespace mat{
    /*
     * 1、规模异常
     */
    class InvalidSizeException: public exception{
    private:
        long long wrong_size;   //错误的规模
        string file_name;       //错误出现的所在文件
        int line_number;        //错误出现的所在行数
        string function_name;   //错误出现的函数名
        string message;         //错误的完整信息
    public:
        //构造器
        InvalidSizeException(long long wrong_size, string file_name, int line_number, string function_name);

        //程序终止时会打印
        [[nodiscard]]
        const char* what () const noexcept override;

        //用户自己做异常处理的时候可以有选择性的打印自己关注的信息，而不用每次都打印what();
        [[nodiscard]]
        long long getWrongSize() const noexcept;

        [[nodiscard]]
        const string &getFileName() const noexcept;

        [[nodiscard]]
        int getLineNumber() const noexcept;

        [[nodiscard]]
        const string &getFunctionName() const noexcept;
    };
    /*
     * 2、空指针异常
     */
    class NullPointerException: public exception{
    private:
        string parameter_name;  //空指针的形参名字
        string file_name;       //错误出现的所在文件
        int line_number;        //错误出现的所在行数
        string function_name;   //错误出现的函数名
        string message;         //错误的完整信息
    public:
        NullPointerException(string parameter_name, string file_name, int line_number, string function_name);

        //程序终止时会打印
        [[nodiscard]]
        const char* what () const noexcept override;

        //用户自己做异常处理的时候可以有选择性的打印自己关注的信息，而不用每次都打印what();
        [[nodiscard]]
        const string &getParameterName() const noexcept;

        [[nodiscard]]
        const string &getFileName() const noexcept;

        [[nodiscard]]
        int getLineNumber() const noexcept;

        [[nodiscard]]
        const string &getFunctionName() const noexcept;
    };
    /*
     * 3、越界异常
     */
    class IndexOutOfBoundException: public exception
    {
    private:
        long long index, bound; //越界的数字和边界大小
        string file_name;       //错误出现的所在文件
        int line_number;        //错误出现的所在行数
        string function_name;   //错误出现的函数名
        string message;         //错误的完整信息
    public:
        IndexOutOfBoundException(long long index, long long bound,
                                 string file_name, int line_number, string function_name);

        //程序终止时会打印
        [[nodiscard]]
        const char* what () const noexcept override;

        //用户自己做异常处理的时候可以有选择性的打印自己关注的信息，而不用每次都打印what();
        [[nodiscard]]
        long long int getIndex() const noexcept;

        [[nodiscard]]
        long long int getBound() const noexcept;

        [[nodiscard]]
        const string &getFileName() const noexcept;

        [[nodiscard]]
        int getLineNumber() const noexcept;

        [[nodiscard]]
        const string &getFunctionName() const noexcept;
    };
    /*
     * 4、运算异常：处理加减法中矩阵size不相等，乘法中第一个列数和第二个的行数不相等
     */
    class OperationException: public exception
    {
    private:
        string fatal_description;   //错误描述
        string file_name;           //错误出现的所在文件
        int line_number;            //错误出现的所在行数
        string function_name;       //错误出现的函数名
        string message;             //错误的完整信息
    public:
        OperationException(string fatal_description,
                                 string file_name, int line_number, string function_name);

        //程序终止时会打印
        [[nodiscard]]
        const char* what () const noexcept override;

        //用户自己做异常处理的时候可以有选择性的打印自己关注的信息，而不用每次都打印what();
        [[nodiscard]]
        const string & getFatalDescription() const noexcept;

        [[nodiscard]]
        const string &getFileName() const noexcept;

        [[nodiscard]]
        int getLineNumber() const noexcept;

        [[nodiscard]]
        const string &getFunctionName() const noexcept;
    };
    /*
     * 5、迭代器越界异常
     */
    class IteratorOutOfBoundException: public exception
    {
    private:
        string fatal_description;   //错误描述
        string file_name;           //错误出现的所在文件
        int line_number;            //错误出现的所在行数
        string function_name;       //错误出现的函数名
        string message;             //错误的完整信息
    public:
        IteratorOutOfBoundException(string fatal_description,
        string file_name, int line_number, string function_name);

        //程序终止时会打印
        [[nodiscard]]
        const char* what () const noexcept override;

        //用户自己做异常处理的时候可以有选择性的打印自己关注的信息，而不用每次都打印what();
        [[nodiscard]]
        const string & getFatalDescription() const noexcept;

        [[nodiscard]]
        const string &getFileName() const noexcept;

        [[nodiscard]]
        int getLineNumber() const noexcept;

        [[nodiscard]]
        const string &getFunctionName() const noexcept;
    };
    /*
     * 6、非方阵导致的越界异常
     */
    class NotASquareMatrixException: public exception
    {
    private:
        size_t rows;                //矩阵的行数
        size_t cols;                //矩阵的列数
        string file_name;           //错误出现的所在文件
        int line_number;            //错误出现的所在行数
        string function_name;       //错误出现的函数名
        string message;             //错误的完整信息
    public:
        NotASquareMatrixException(size_t rows, size_t cols,
                                    string file_name, int line_number, string function_name);

        //程序终止时会打印
        [[nodiscard]]
        const char* what () const noexcept override;

        //用户自己做异常处理的时候可以有选择性的打印自己关注的信息，而不用每次都打印what();
        [[nodiscard]]
        const size_t & getRows() const noexcept;

        [[nodiscard]]
        const size_t & getCols() const noexcept;

        [[nodiscard]]
        const string &getFileName() const noexcept;

        [[nodiscard]]
        int getLineNumber() const noexcept;

        [[nodiscard]]
        const string &getFunctionName() const noexcept;
    };
    /*
     * 7、范围不合法
     */
    class InvalidRangeException: public exception
    {
    private:
        size_t lb;                //范围的下界
        size_t ub;                //范围的上界
        string file_name;           //错误出现的所在文件
        int line_number;            //错误出现的所在行数
        string function_name;       //错误出现的函数名
        string message;             //错误的完整信息
    public:
        InvalidRangeException(size_t lb, size_t ub,
                              string file_name, int line_number, string function_name);

        //程序终止时会打印
        [[nodiscard]]
        const char* what () const noexcept override;

        //用户自己做异常处理的时候可以有选择性的打印自己关注的信息，而不用每次都打印what();
        [[nodiscard]]
        const size_t & getUb() const noexcept;

        [[nodiscard]]
        const size_t & getLb() const noexcept;

        [[nodiscard]]
        const string &getFileName() const noexcept;

        [[nodiscard]]
        int getLineNumber() const noexcept;

        [[nodiscard]]
        const string &getFunctionName() const noexcept;
    };

}



#endif //CPPPROJECTS_MATEXCEPTION_HPP
