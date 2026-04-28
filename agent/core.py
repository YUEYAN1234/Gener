"""
Agent 核心引擎 - Core Engine
对话管理 + 工具调度循环 + 流式逐字输出
"""

import json
import sys
import time
from typing import Optional

from agent.deepseek_client import DeepSeekClient
from agent.system_prompt import SYSTEM_PROMPT
from agent.tool_registry import get_all_tools, execute_tool, get_tool_names
from tools.biosafety_checker import BiosafetyViolation


# 逐字打印速度 (秒/字符)
CHAR_DELAY = 0.015
FAST_CHAR_DELAY = 0.005


def _type_print(text: str, delay: float = CHAR_DELAY, end: str = "\n"):
    """逐字打印效果"""
    for ch in text:
        sys.stdout.write(ch)
        sys.stdout.flush()
        time.sleep(delay)
    if end:
        sys.stdout.write(end)
        sys.stdout.flush()


def _truncate_args(args: dict, max_len: int = 80) -> str:
    """截断参数显示"""
    s = json.dumps(args, ensure_ascii=False)
    if len(s) > max_len:
        return s[:max_len] + "..."
    return s


class DNADesignAgent:
    """大片段 DNA 进化设计与物化分析智能体"""

    def __init__(self):
        self.client = DeepSeekClient()
        self.tools = get_all_tools()
        self.messages = [
            {"role": "system", "content": SYSTEM_PROMPT},
        ]
        self.max_tool_rounds = 20  # 提升到 20 轮，避免复杂任务中断

    def chat(self, user_input: str) -> str:
        """
        处理用户输入，工具调用阶段显示思考过程，最终回复逐字流式输出。

        流程: 用户消息 -> DeepSeek -> [tool_calls -> 执行 -> 回传] x N -> 流式最终回复
        """
        self.messages.append({"role": "user", "content": user_input})

        for round_idx in range(self.max_tool_rounds):
            # 分隔线
            print(f"\n{'- ' * 30}")
            if round_idx == 0:
                _type_print("[思考中] 正在理解需求并规划分析路径...", FAST_CHAR_DELAY)
            else:
                _type_print(f"[推理] 第 {round_idx + 1} 轮工具调用...", FAST_CHAR_DELAY)

            try:
                response = self.client.chat(
                    messages=self.messages,
                    tools=self.tools,
                )
            except Exception as e:
                error_msg = f"[错误] API 调用失败: {e}"
                print(error_msg)
                return error_msg

            choice = response.choices[0]
            message = choice.message

            # 如果模型产生了思考内容 (reasoning_content)，展示出来
            reasoning = getattr(message, "reasoning_content", None)
            if reasoning:
                print(f"\n  [内部推理]")
                for line in reasoning.strip().split("\n"):
                    _type_print(f"  | {line}", FAST_CHAR_DELAY)

            # 检查是否有工具调用
            if message.tool_calls:
                self.messages.append(message.model_dump())

                # 如果有文本内容伴随工具调用，先展示
                if message.content:
                    print(f"\n  [分析思路]")
                    for line in message.content.strip().split("\n"):
                        _type_print(f"  > {line}", FAST_CHAR_DELAY)

                for tool_call in message.tool_calls:
                    func_name = tool_call.function.name
                    try:
                        arguments = json.loads(tool_call.function.arguments)
                    except json.JSONDecodeError:
                        arguments = {}

                    _type_print(f"  [调用] {func_name}({_truncate_args(arguments)})", FAST_CHAR_DELAY)

                    try:
                        result = execute_tool(func_name, arguments)
                    except BiosafetyViolation as bsv:
                        alert = (
                            f"\n{'!' * 60}\n"
                            f"  [生物安全警报] 检测到管制病原体序列!\n"
                            f"  病原体: {bsv.pathogen}\n"
                            f"  基因: {bsv.gene}\n"
                            f"  风险等级: Risk Group {bsv.risk_group}\n"
                            f"  相似度: {bsv.identity:.1%}\n"
                            f"\n  任务已终止。根据生物安全法规，禁止合成此类序列。\n"
                            f"{'!' * 60}"
                        )
                        print(alert)
                        self.messages.append({
                            "role": "tool",
                            "tool_call_id": tool_call.id,
                            "content": json.dumps({"BIOSAFETY_VIOLATION": str(bsv)}),
                        })
                        return alert
                    except Exception as e:
                        result = {"error": f"{type(e).__name__}: {str(e)}"}

                    result_str = json.dumps(result, ensure_ascii=False, default=str)
                    if len(result_str) > 30000:
                        result_str = result_str[:30000] + '... [TRUNCATED]'

                    _type_print(f"  [完成] {func_name} -> {len(result_str)} chars", FAST_CHAR_DELAY)

                    self.messages.append({
                        "role": "tool",
                        "tool_call_id": tool_call.id,
                        "content": result_str,
                    })

                continue

            else:
                # 没有工具调用 -> 最终回复，使用流式逐字输出
                print(f"\n{'- ' * 30}")
                _type_print("[生成回复]\n", FAST_CHAR_DELAY)
                reply = message.content or ""
                self.messages.append({"role": "assistant", "content": reply})
                return reply

        return "[警告] 达到最大工具调用轮数限制(20)。请尝试简化请求。"

    def chat_stream(self, user_input: str) -> str:
        """
        流式处理：工具调用阶段非流式，最终回复流式逐字输出。
        """
        self.messages.append({"role": "user", "content": user_input})

        for round_idx in range(self.max_tool_rounds):
            print(f"\n{'- ' * 30}")
            if round_idx == 0:
                _type_print("[思考中] 正在理解需求并规划分析路径...", FAST_CHAR_DELAY)
            else:
                _type_print(f"[推理] 第 {round_idx + 1} 轮工具调用...", FAST_CHAR_DELAY)

            try:
                # 工具调用阶段用非流式
                response = self.client.chat(
                    messages=self.messages,
                    tools=self.tools,
                    stream=False,
                )
            except Exception as e:
                return f"[错误] API 调用失败: {e}"

            choice = response.choices[0]
            message = choice.message

            # 展示推理过程
            reasoning = getattr(message, "reasoning_content", None)
            if reasoning:
                print(f"\n  [内部推理]")
                for line in reasoning.strip().split("\n"):
                    _type_print(f"  | {line}", FAST_CHAR_DELAY)

            if message.tool_calls:
                self.messages.append(message.model_dump())

                if message.content:
                    print(f"\n  [分析思路]")
                    for line in message.content.strip().split("\n"):
                        _type_print(f"  > {line}", FAST_CHAR_DELAY)

                for tool_call in message.tool_calls:
                    func_name = tool_call.function.name
                    try:
                        arguments = json.loads(tool_call.function.arguments)
                    except json.JSONDecodeError:
                        arguments = {}

                    _type_print(f"  [调用] {func_name}({_truncate_args(arguments)})", FAST_CHAR_DELAY)

                    try:
                        result = execute_tool(func_name, arguments)
                    except BiosafetyViolation as bsv:
                        alert = f"\n[生物安全警报] {bsv}\n任务已终止。"
                        self.messages.append({
                            "role": "tool",
                            "tool_call_id": tool_call.id,
                            "content": json.dumps({"BIOSAFETY_VIOLATION": str(bsv)}),
                        })
                        return alert
                    except Exception as e:
                        result = {"error": str(e)}

                    result_str = json.dumps(result, ensure_ascii=False, default=str)
                    if len(result_str) > 30000:
                        result_str = result_str[:30000] + '... [TRUNCATED]'

                    _type_print(f"  [完成] {func_name} -> {len(result_str)} chars", FAST_CHAR_DELAY)

                    self.messages.append({
                        "role": "tool",
                        "tool_call_id": tool_call.id,
                        "content": result_str,
                    })
                continue

            else:
                # 最终回复 -> 真正的流式逐字输出
                print(f"\n{'- ' * 30}")
                _type_print("[生成回复]\n", FAST_CHAR_DELAY)

                try:
                    stream = self.client.chat(
                        messages=self.messages,
                        tools=None,  # 最终回复不再带工具
                        stream=True,
                    )
                    full_reply = ""
                    for chunk in stream:
                        delta = chunk.choices[0].delta
                        if delta.content:
                            for ch in delta.content:
                                sys.stdout.write(ch)
                                sys.stdout.flush()
                                time.sleep(CHAR_DELAY)
                            full_reply += delta.content
                    print()  # 换行
                    self.messages.append({"role": "assistant", "content": full_reply})
                    return full_reply
                except Exception as e:
                    # fallback: 非流式
                    reply = message.content or ""
                    self.messages.append({"role": "assistant", "content": reply})
                    return reply

        return "[警告] 达到最大工具调用轮数限制(20)。请尝试简化请求。"

    def reset(self):
        """重置对话历史"""
        self.messages = [
            {"role": "system", "content": SYSTEM_PROMPT},
        ]
        print("[系统] 对话已重置。")

    def get_tool_list(self) -> list[str]:
        """获取所有可用工具名称"""
        return get_tool_names()
