"""
DeepSeek API 客户端 - DeepSeek Client
基于 OpenAI SDK 的 DeepSeek 封装，支持 function calling + 流式输出
"""

import os
import json
from openai import OpenAI
from dotenv import load_dotenv

load_dotenv()


class DeepSeekClient:
    def __init__(self):
        self.api_key = os.getenv("DEEPSEEK_API_KEY")
        self.base_url = os.getenv("DEEPSEEK_BASE_URL", "https://api.deepseek.com")
        self.model = os.getenv("DEEPSEEK_MODEL", "deepseek-v4-pro")

        if not self.api_key:
            raise ValueError("DEEPSEEK_API_KEY not found in environment variables")

        self.client = OpenAI(
            api_key=self.api_key,
            base_url=self.base_url,
        )

    def chat(
        self,
        messages: list,
        tools: list = None,
        temperature: float = 0.3,
        top_p: float = 0.9,
        max_tokens: int = 8192,
        stream: bool = False,
    ):
        """
        发送对话请求。

        Args:
            messages: 消息列表 [{role, content}]
            tools: function calling 工具定义
            temperature: 温度 (0.3 for scientific rigor)
            top_p: Top-p 采样
            max_tokens: 最大输出 token
            stream: 是否流式

        Returns:
            ChatCompletion 或 Stream
        """
        kwargs = {
            "model": self.model,
            "messages": messages,
            "temperature": temperature,
            "top_p": top_p,
            "max_tokens": max_tokens,
            "stream": stream,
        }

        if tools:
            kwargs["tools"] = tools
            kwargs["tool_choice"] = "auto"

        response = self.client.chat.completions.create(**kwargs)
        return response

    def chat_stream(self, messages, tools=None, temperature=0.3, top_p=0.9, max_tokens=8192):
        """流式对话，返回生成器"""
        return self.chat(messages, tools, temperature, top_p, max_tokens, stream=True)
