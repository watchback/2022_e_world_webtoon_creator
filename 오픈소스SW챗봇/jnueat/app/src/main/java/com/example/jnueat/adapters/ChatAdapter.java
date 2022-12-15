package com.example.jnueat.adapters;

import android.app.Activity;
import android.speech.tts.TextToSpeech;
import android.view.LayoutInflater;
import android.view.View;
import android.view.ViewGroup;
import android.widget.TextView;

import androidx.annotation.NonNull;
import androidx.recyclerview.widget.RecyclerView;

import com.example.jnueat.R;
import com.example.jnueat.models.Message;

import java.util.List;

public class ChatAdapter extends RecyclerView.Adapter<ChatAdapter.MyViewHolder> {

  private List<Message> messageList;
  private Activity activity;

  //tts
  private TextToSpeech mTTS;

  public ChatAdapter(List<Message> messageList, Activity activity) {
    this.messageList = messageList;
    this.activity = activity;
  }

  @NonNull @Override
  public MyViewHolder onCreateViewHolder(@NonNull ViewGroup parent, int viewType) {
    View view = LayoutInflater.from(activity).inflate(R.layout.adapter_message_one, parent, false);

    return new MyViewHolder(view);
  }

  @Override public void onBindViewHolder(@NonNull final MyViewHolder holder, int position) {
    final String message = messageList.get(position).getMessage();
    boolean isReceived = messageList.get(position).getIsReceived();

    //드래그
    holder.messageReceive.setTextIsSelectable(false);
    holder.messageReceive.measure(-1, -1);
    holder.messageReceive.setTextIsSelectable(true);

     if(isReceived){
       holder.messageReceive.setVisibility(View.VISIBLE);
       holder.messageSend.setVisibility(View.GONE);
       holder.messageReceive.setText(message);
     }else {
       holder.messageSend.setVisibility(View.VISIBLE);
       holder.messageReceive.setVisibility(View.GONE);
       holder.messageSend.setText(message);
     }
  }

  @Override public int getItemCount() {
    return messageList.size();
  }

  static class MyViewHolder extends RecyclerView.ViewHolder{

    TextView messageSend;
    TextView messageReceive;

    MyViewHolder(@NonNull View itemView) {
      super(itemView);
      messageSend = itemView.findViewById(R.id.message_send);
      messageReceive = itemView.findViewById(R.id.message_receive);
    }
  }

}
