package com.example.jnueat;

import androidx.appcompat.app.AppCompatActivity;

import android.content.Intent;
import android.os.Bundle;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.Toast;

import com.android.volley.RequestQueue;
import com.android.volley.Response;
import com.android.volley.toolbox.Volley;

import org.json.JSONException;
import org.json.JSONObject;

/**
 * 로그인 화면과 연계하여 로그인을 처리하는 클래스이다.
 * @author 김한진
 */

public class loginActivity extends AppCompatActivity {
    private EditText ett_id, ett_pass;
    private Button btn_login, btn_register2;


    /**
     * 로그인화면의 버튼들과 아이디 패스워드 입력하는 칸을 연결 해주는 액티비티 메소드
     * @param savedInstanceState
     */
    @Override
    protected void onCreate(Bundle savedInstanceState) { //액티비티 시작시 처음으로 실행되는 생명주기
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_login);

        ett_id = findViewById(R.id.ett_id); // 아이디 빈칸 연결하는 객체
        ett_pass = findViewById(R.id.ett_pass); // 패스워드 빈칸 연결하는 객체
        btn_login = findViewById(R.id.btn_login); // 로그인 버튼 연결하는 객체

        //로그인 버튼 클릭 시 수행
        btn_register2 = findViewById(R.id.btn_register2); // 회원가입 버튼 연결하는 객체
        btn_register2.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) { // EditText에 현재 입력되어 있는 값을 get(가져온다)
                Intent intent = new Intent(loginActivity.this,registerActivity.class);
                startActivity(intent);
            }
        }); //회원가입 버튼 누르면 발생하는 액티비티 설정

        btn_login.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                String userID = ett_id.getText().toString();
                String userPass = ett_pass.getText().toString();

                Response.Listener<String> responseListener = new Response.Listener<String>() {
                    @Override
                    public void onResponse(String response) {
                        try {
                            JSONObject jsonObject = new JSONObject(response);
                            boolean success = jsonObject.getBoolean("success");
                            if (success) { //로그인 성공
                                String userID = jsonObject.getString("userID");
                                String userPass = jsonObject.getString("userPassword");
                                Toast.makeText(getApplicationContext(), "로그인에 성공", Toast.LENGTH_SHORT).show();
                                Intent intent = new Intent(loginActivity.this, MainActivity.class);
                                intent.putExtra("userID", userID);
                                intent.putExtra("userPass", userPass);
                                startActivity(intent);
                            } else { //로그인에 실패
                                Toast.makeText(getApplicationContext(), "로그인에 실패", Toast.LENGTH_SHORT).show();
                                return;
                            }
                        } catch (JSONException e) {
                            e.printStackTrace();
                        }


                    }
                };
                // 서버로 Volley를 이용해서 요청을 함
                LoginRequest loginRequest = new LoginRequest(userID, userPass, responseListener);
                RequestQueue queue = Volley.newRequestQueue(loginActivity.this);
                queue.add(loginRequest);


            }
        });// 로그인 버튼 누르면 발생하는 액티비티 설정

    }
}