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
 * 회원가입 화면과 연계하여 회원가입을 처리하는 클래스이다.
 * @author 김한진
 */
public class registerActivity extends AppCompatActivity {

    private EditText et_id, et_pass, et_name, et_age;
    private Button btn_register;

    /**
     * 회원가입 화면의 버튼들과 아이디 패스워드 및 이름,나이 입력하는 칸을 연결 해주는 액티비티 메소드
     * @param savedInstanceState
     */
    @Override
    protected void onCreate(Bundle savedInstanceState) { //액티비티 시작시 처음으로 실행되는 생명주기
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_register);

        // 아이디 값 찾아주기
        et_id = findViewById(R.id.et_id); // 아이디 빈칸 연결하는 객체
        et_pass = findViewById(R.id.et_pass);// 비밀번호 빈칸 연결하는 객체
        et_name = findViewById(R.id.et_name);// 이름 빈칸 연결하는 객체
        et_age = findViewById(R.id.et_age);// 나이 빈칸 연결하는 객체

        //회원가입 버튼 클릭 시 수행
        btn_register = findViewById(R.id.btn_register);
        btn_register.setOnClickListener(new View.OnClickListener() {

            /**
             * EditText에 현재 입력되어 있는 값을 get(가져온다)
             * @param view
             */
            @Override
            public void onClick(View view) {
                String userID = et_id.getText().toString();
                String userPass = et_pass.getText().toString();
                String userName = et_name.getText().toString();
                int userAge = Integer.parseInt(et_age.getText().toString());

                Response.Listener<String> responseListener = new Response.Listener<String>() {


                    /**
                     * 회원가입 요청하는 메소드
                     * @param response
                     */
                    @Override
                    public void onResponse(String response) {
                        try {
                            JSONObject jsonObject = new JSONObject(response);
                            boolean success = jsonObject.getBoolean("success");
                            if (success) {
                                Toast.makeText(getApplicationContext(), "회원 등록에 성공", Toast.LENGTH_SHORT).show();
                                Intent intent = new Intent(registerActivity.this, loginActivity.class);
                                startActivity(intent);
                            } else { //회원등록에 실패
                                Toast.makeText(getApplicationContext(), "회원 등록에 실패", Toast.LENGTH_SHORT).show();
                                return;
                            }
                        } catch (JSONException e) {
                            e.printStackTrace();
                        }

                    }
                };
                // 서버로 Volley를 이용해서 요청을 함
                RegisterRequest registerRequest = new RegisterRequest(userID,userPass,userName,userAge,responseListener);
                RequestQueue queue = Volley.newRequestQueue(registerActivity.this);
                queue.add(registerRequest);

            }
        });
    }
}
