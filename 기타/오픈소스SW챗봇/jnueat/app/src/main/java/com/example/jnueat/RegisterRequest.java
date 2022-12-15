package com.example.jnueat;

import com.android.volley.AuthFailureError;
import com.android.volley.Response;
import com.android.volley.toolbox.StringRequest;

import java.util.HashMap;
import java.util.Map;

/**
 * 서버에 회원가입 요청을 하는 클래스이다.
 * @author 김한진
 */
public class RegisterRequest extends StringRequest {

    //서버 URL 설정 (PHP 파일 연동)
    final static private String URL = "http://oxox6300.dothome.co.kr/Register.php";
    private Map<String, String> map;

    /**
     * 유저 아이디,패스워드,이름,나이를 데이터 베이스에 키 값을 비교해서 삽입한다.
     * @param userID
     * @param userPassword
     * @param userName
     * @param userAge
     * @param listener
     */
    public RegisterRequest(String userID, String userPassword, String userName, int userAge, Response.Listener<String> listener) {
        super(Method.POST, URL, listener, null);

        map = new HashMap<>();
        map.put("userID",userID);
        map.put("userPassword", userPassword);
        map.put("userName", userName);
        map.put("userAge", userAge +"");
    }

    @Override
    protected Map<String, String> getParams() throws AuthFailureError {
        return map;
    }
}
