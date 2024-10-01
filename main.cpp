#include "bswabe.h"
#include <pbc/pbc.h>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <algorithm>

using namespace std;
const int ELEMENT_SIZE = 2048;        // 假设曲线元素的字节大小
// 将字符串转换为element
void string_to_element(element_t& element, const string& str) {
    int len = str.size();
    // 确保缓冲区包含固定大小并为字符串文字和填充可能更长的元素而设计
    vector<unsigned char> buffer(ELEMENT_SIZE, 0);
    memcpy(buffer.data(), str.c_str(), len);
    element_from_bytes(element, buffer.data());
}
// 将element转换为字符串
void element_to_string(element_t& element, string& str) {
    int length = element_length_in_bytes(element);
    vector<unsigned char> buffer(length);
    element_to_bytes(buffer.data(), element);
    // 从缓冲区中提取原始字符串内容
    auto pos = std::find(buffer.begin(), buffer.end(), 0);
    str.assign(reinterpret_cast<char*>(buffer.data()), pos != buffer.end() ? pos - buffer.begin() : length);
}

int main() {
    cout<<endl<<"********************************************step1:系统初始化阶段***************************************************"<<endl;
    BswabePub pub;             // 定义系统公钥
    BswabeMsk msk;             // 定义系统主私钥
    bswabe_setup(&pub, &msk);  // 初始化系统，生成公钥和主私钥

    cout<<endl<<"********************************************step2:密钥生成阶段***************************************************"<<endl;
    std::vector<std::string> user_attrs = {"attr25", "attr6", "attr7s","attr1","attr7","attr83" };   // 定义用户属性集
    // std::vector<std::string> user_attrs = {"attr4","attr2","attr3", "attr6"};   // 定义用户属性集
    BswabePrv* prv = bswabe_keygen(&pub, &msk, user_attrs);   // 生成用户私钥

    cout<<endl<<"********************************************step3:加密阶段***************************************************"<<endl;
    std::string policy = "attr1 attr2 attr3 attr4 1of2 2of2 1of2 attr5 attr6 attr7 attr8 1of2 2of2 1of2 2of2";   // 定义策略树
    // std::string policy = "attr1 attr2 attr3 attr4 3of4";   // 定义策略树
    
    // 当前曲线最大加密长度为128 bytes，现在输入是140 bytes，但是后面解密只能恢复前面128 bytes
    string M_string = "cc11111111112222222222333333333344444444445555555555666666666677777777778aabbbbbbbbbbccccccccccddddddddddeeeeeeeeee33";
    element_t M;                        // element_t明文
    element_init_GT(M, pub.pairing);    // 初始化M
    string_to_element(M, M_string);     // 将string转成element
    std::cout << "Plaintext M_string: " << M_string << ", size = "<<M_string.size()<<std::endl;
    element_printf("Plaintext M: %B\n", M);
    BswabeCph* cph = bswabe_enc(&pub, policy, M);   // 使用policy对M进行加密
    if (cph == nullptr) {
        std::cerr << "Encrypt failed!" << std::endl;
        return 1;
    }
    element_printf("cs (Encrypted message): %B\n", cph->cs);    // 打印加密后的明文 \tilde{C}
    std::cout << "Access policy tree:" << std::endl;
    print_policy_tree(cph->p, 0);                       // 打印访问策略树的内容

    cout<<endl<<"********************************************step4:解密阶段***************************************************"<<endl;
    BswabeElementBoolean* result = bswabe_dec(&pub, prv, cph);   // 使用用户私钥prv对密文cph进行解密
    if (result->b) {
        string output;
        element_to_string(result->e, output);       // 将element转成string
        cout << "Decrypt success! Decrypted plaintext is: " << output <<",size= "<<output.size()<<endl;   // 解密成功
    } else {
        std::cout << "Decrypt failed!" << std::endl;
    }

    // 释放资源
    free_bswabe_prv(prv);
    free_bswabe_cph(cph);
    element_clear(result->e);
    delete result;
    element_clear(msk.beta);
    element_clear(msk.g_alpha);
    element_clear(pub.g);
    element_clear(pub.h);
    element_clear(pub.f);
    element_clear(pub.e_gg_alpha);
    pairing_clear(pub.pairing);

    return 0;
}
