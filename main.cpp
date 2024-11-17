#include <iostream>
#include <memory>
#include "matrix.hpp"
#include "MATexception.hpp"
using namespace std;


int main(){
//    //1.3.2 verification
//    //构造器1
//    mat::matI32 m1(2, 2, 3);
//    cout<<m1<<endl;
//    //构造器2
//    int arr[4] = {1,2,3,4};
//    mat::matI32 m2(2, 2, 1, arr);
//    cout<<m2<<endl;
//    //构造器3: 双通道的初始化列表
//    mat::matI32 m3(2, 2, 2,
//                   {1, 2,   3, 4,
//                    5, 6,   7, 8});
//    cout<<m3<<endl;
//    //构造器4
//    mat::matI32 m4(2, 2, 1, 5);
//    cout<<m4<<endl;
//    //构造器5
//    mat::matI32 m5(3, 1, mat::mat_type::identity_mat);
//    cout<<m5<<endl;
//    //构造器6
//    mat::matI32 m6(m5);
//    cout<<m6<<endl;

//    //2.2.3 verification
//    mat::matI32 m1(2, 2, 2,
//                   {1, 2,   3, 4,
//                    5, 6,   7, 8});
//    cout<<"print m1:"<<endl;
//    cout<<m1<<endl;
//    cout<<"modify m1:"<<endl;
//    for(auto it = m1.getBegin(); it != m1.getEnd(); it++)
//    {
//        *it = 100;
//    }
//    cout<<"print m1 through iterator:"<<endl;
//    for(auto it = m1.getBegin(); it != m1.getEnd(); it++)
//    {
//        cout<<*it<<endl;
//    }

//    //3.4 verification
//    mat::matI32 m1(5, 5, 3,
//                   {1,1,1,  1,2,1,  1,3,1,  1,4,1,  1,5,1,
//                    1,2,1,  1,3,1,  1,4,1,  1,5,1,  1,6,1,
//                    1,3,1,  1,4,1,  1,5,1,  1,6,1,  1,7,1,
//                    1,4,1,  1,5,1,  1,6,1,  1,7,1,  1,8,1,
//                    1,5,1,  1,6,1,  1,7,1,  1,8,1,  1,9,1});
//    cout<<"print m1: "<<endl;
//    cout<<m1<<endl;
//    cout<<"m1 add a new channel: "<<endl;
//    m1.add_channel({2, 2, 2, 2, 2,
//                    2, 2, 2, 2, 2,
//                    2, 2, 2, 2, 2,
//                    2, 2, 2, 2, 2,
//                    2, 2, 2, 2, 2});
//    cout<<m1<<endl;
//    cout<<"m1 erase the index 2 channel(the 3rd channel): "<<endl;
//    m1.erase_channel(2);
//    cout<<m1<<endl;
//    //m2 是 m1 第3行，第3列，第2层截取出来的子矩阵
//    auto m2 = mat::matrix<int>(m1, mat::mat_range(2, 3), mat::mat_range(2, 3), mat::mat_range(1, 2));
//    cout<<"print m2:"<<endl;
//    cout<<m2<<endl;
//    cout<<"Enlarge m2 2 units in row dimension, 2 units in col dimension, 1 units in channel dimension:"<<endl;
//    m2.reshape(2, 2, 1);
//    cout<<m2<<endl;
//    cout<<"Shrink m2 1 unit in row dimension, 1 unit in col dimension, 1 units in channel dimension:"<<endl;
//    m2.reshape(-1, -1, -1);
//    cout<<m2<<endl;
//    cout<<"Move m2 1 units to the right, and 1 unit downward: "<<endl;
//    m2.move(1, 1, 0);
//    cout<<m2<<endl;
//    cout<<"Judge whether m2 is a sub-matrix of m1: "<<endl;
//    cout<<(m2.isSubMatrix(m1)?"is a sub-matrix":"not a sub-matrix")<<endl;

//    //4.1.6 verification
//    mat::matI32 m1(3, 3, 2,
//                   {0,1, 0,1, 0,1,
//                    0,1, 0,1, 0,1,
//                    0,1, 0,1, 0,1});
//
//    cout<<"assign to m2"<<endl;
//    mat::matI32 m2;
//    m2 = m1;
//    cout<<m2<<endl;
//
//    cout<<"print 2nd channel: "<<endl;
//    for(int i = 0; i<m1.getRows(); i++)
//    {
//        for(int j = 0; j<m1.getCols(); j++)
//            cout<<m1[1][i * m1.getStep() + j]<<" ";
//        cout<<endl;
//    }
//    cout<<endl;
//
//    cout<<"change channel 2, row 2, col 3 to 5"<<endl;
//    m1(1, 1, 2) = 5;
//    cout<<m1<<endl;
//
//    cout<<"print the fields:"<<endl;
//    cout<<"rows: "<<m2.getRows()<<" cols: "<<m2.getCols()<<" channels: "<<m2.getChannels()
//    <<" shift: "<<m2.getShift()<<" step: "<<m2.getStep()<<" step_row: "<<m2.getStepRows()<<endl<<endl;
//
//    cout<<"set the whole channel1 to 3"<<endl;
//    m1.set_channel(0, 3);
//    cout<<m1<<endl;
//
//    cout<<"set channel1, row2 to 1"<<endl;
//    m1.set_row(0, 1, 1);
//    cout<<m1<<endl;
//
//    cout<<"set channel2, col2 to 10"<<endl;
//    m1.set_col(1, 1, 10);
//    cout<<m1<<endl;

//    //4.2.1.2 verification
//    mat::matI32 m1(3, 3, 1,{0, 0, 0,
//                            0, 0, 0,
//                            0, 0, 0});
//    //截取第二行
//    cout<<"the second row: "<<endl;
//    mat::matI32 m2(m1, mat::mat_range(1, 2), mat::mat_range(0, 3), mat::mat_range(0, 1));
//    cout<<m2<<endl;
//
//    //m3 = m2 + 3;
//    cout<<"m3 = m2 + 3: "<<endl;
//    mat::matI32 m3 = m2 + 3;
//    cout<<m3<<endl;
//
//    //m2自己加1并且打印m1
//    cout<<"m2 += 1: "<<endl;
//    m2 += 1;
//    cout<<m1<<endl;
//
//    //m2自己加上m3并且打印m2
//    cout<<"m2 += m3: "<<endl;
//    m2 += m3;
//    cout<<m2<<endl;

//    //4.2.2.2 verification
//    mat::matI32 m1(3, 3, 1,{1, 1, 1,
//                            1, 1, 1,
//                            1, 1, 1});
//    mat::matI32 m2 = 3 - m1;
//    cout<<m2<<endl;

//    //4.2.3.2 verification
//    mat::matI32 m1(3, 3, 1,{1, 1, 1,
//                            1, 1, 1,
//                            1, 1, 1});
//    //截取两个子矩阵进行相乘
//    auto m2 = mat::matI32(m1, mat::mat_range(0, 2), mat::mat_range(0, 2), mat::mat_range(0, 1)) *
//              mat::matI32(m1, mat::mat_range(1, 3), mat::mat_range(1, 3), mat::mat_range(0, 1));
//    cout<<m2<<endl;

//    //4.3.2 verification
//    mat::matI32 m1(3, 3, 1, {1, 0, 1,
//                             0, 1, 0,
//                             1, 0, 1});
//
//    mat::matI32 m2(3, 3, 1, {0, 1, 0,
//                             1, 0, 1,
//                             0, 1, 0});
//
//    cout<<"^ operation:"<<endl;
//    cout<<(m1 ^ m2)<<endl;
//    cout<<"& operation:"<<endl;
//    cout<<(m1 & m2)<<endl;
//    cout<<"| operation:"<<endl;
//    cout<<(m1 | m2)<<endl;

//    //4.4.2 verification
//    mat::matI32 m1(2, 4, 1, {0, 0, 1, 1,
//                             0, 0, 1, 1});
//    mat::matI32 m2(2, 4, 1, {0, 0, 2, 2,
//                             0, 0, 2, 2});
//    cout<<"== operation:"<<endl;
//    cout<<(mat::matI32(m1, mat::mat_range(0, 2), mat::mat_range(0, 2), mat::mat_range(0, 1))
//          ==
//          mat::matI32(m2, mat::mat_range(0, 2), mat::mat_range(0, 2), mat::mat_range(0, 1)))<<endl;
//    cout<<"< operation:"<<endl;
//    cout<<(mat::matI32(m1, mat::mat_range(0, 2), mat::mat_range(2, 4), mat::mat_range(0, 1))
//           <
//           mat::matI32(m2, mat::mat_range(0, 2), mat::mat_range(2, 4), mat::mat_range(0, 1)))<<endl;

//    //4.5.2 verification
//    mat::matI32 m1(2, 2, 1, {0, 0,
//                             0, 0});
//    mat::matI32 m2(2, 2, 1, {0, 0,
//                             0, 0});
//    cout<<"print ++m1:"<<endl;
//    cout<<(++m1)<<endl;
//    cout<<"print m2++:"<<endl;
//    cout<<(m2++)<<endl;

//    //4.6.2 verification
//    mat::matI32 m1(2, 2, 1, {7, 7,
//                             7, 7});
//    cout<<(m1<<2)<<endl;
//    cout<<(m1>>2)<<endl;

//    //4.7.2 verification
//    mat::matI32 m1(3, 3, 2);
//    cin>>m1;
//    cout<<m1;
//    m1.set_IO_mode(mat::io_mode::row_by_row_col_by_col);
//    cout<<m1;

    //5 verification
//    mat::matI32 m1(2, 2, 1, {1, 1,
//                             1, 1});
//    mat::matI32 m2(1, 1, 1, {5});
//    cout<<(m1 * m2)
}


