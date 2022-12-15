package com.example.jnueat;

import com.android.volley.AuthFailureError;
import com.android.volley.Response;
import com.android.volley.toolbox.StringRequest;

import java.util.HashMap;
import java.util.Map;

/**
 * 서버에 로그인 요청을 하는 클래스이다.
 * @author 김한진
 */
public class LoginRequest extends StringRequest {

    //서버 URL 설정 (PHP 파일 연동)
    final static private String URL = "http://oxox6300.dothome.co.kr/Login.php";
    private Map<String, String> map;

    /**
     * 유저 아이디, 패스워드의 키 값을 데이터 베이스와 비교하는 메소드
     * @param userID
     * @param userPassword
     * @param listener
     */

    public LoginRequest(String userID, String userPassword, Response.Listener<String> listener) {
        super(Method.POST, URL, listener, null);

        map = new HashMap<>();
        map.put("userID",userID);
        map.put("userPassword", userPassword);
    }

    @Override
    protected Map<String, String> getParams() throws AuthFailureError {
        return map;
    }
}
