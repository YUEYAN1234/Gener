"""
大片段 DNA 进化设计与物化分析智能体 Gener - CLI 入口
"""

import sys
import os
import io
import time

# Windows 下强制 UTF-8 输出
if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

# 确保项目根目录在 path 中
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from agent.core import DNADesignAgent, _type_print, FAST_CHAR_DELAY


BANNER = r"""
  ____                          
 / ___| ___ _ __   ___ _ __     
| |  _ / _ \ '_ \ / _ \ '__|   
| |_| |  __/ | | |  __/ |      
 \____|\___|_| |_|\___|_|       

 Large-Fragment DNA Evolutionary Design
 & Physicochemical Analysis Agent
"""

DIVIDER_THICK = "=" * 56
DIVIDER_THIN = "-" * 56

INFO_BLOCK = f"""
{DIVIDER_THICK}
  Gener - 大片段 DNA 进化设计与物化分析智能体
{DIVIDER_THICK}
  核心引擎    DeepSeek API (Function Calling)
  工具数      {{tool_count}} 个生物学分析工具
  能力域      序列分析 | 热力学计算 | 结构预测
              loxP/SCRaMbLE | 端粒/着丝粒设计
              Gibson/Golden Gate 拆分 | 生物安全
{DIVIDER_THIN}
  命令:
    /tools    列出所有可用工具
    /reset    重置对话上下文
    /help     显示使用帮助
    /quit     退出程序
{DIVIDER_THICK}
"""

HELP_TEXT = f"""
{DIVIDER_THIN}
  使用指南
{DIVIDER_THIN}

  1. 直接输入自然语言描述你的设计需求
  2. 可以粘贴 DNA 序列进行分析
  3. 支持多行输入，连按两次 Enter 发送
  4. 支持多轮对话，Agent 会记住上下文

  示例输入:
    "分析这段序列的合成难度: ATCGATCG..."
    "设计一段 80kb 的合成酵母染色体臂，包含端粒和 loxP 位点"
    "帮我把这段 20kb 的序列拆分成 Gibson Assembly 片段"
    "检查这段序列中是否有 G-quadruplex 结构"

  工具涵盖:
    序列验证, 热力学计算, G4/Z-DNA结构分析, 
    合成难度评分, loxP设计, SCRaMbLE模拟,
    端粒/着丝粒设计, Gibson/Golden Gate片段拆分,
    生物安全筛查, GenBank/FASTA格式导出

{DIVIDER_THIN}
"""


def main():
    # 打印 banner（逐行快速显示）
    for line in BANNER.strip().split("\n"):
        print(f"  {line}")
        time.sleep(0.05)
    print()

    # 初始化 Agent
    try:
        agent = DNADesignAgent()
        info = INFO_BLOCK.format(tool_count=len(agent.tools))
        print(info)
        _type_print(f"  [就绪] 模型: {agent.client.model} | 工具: {len(agent.tools)} 个已加载", FAST_CHAR_DELAY)
        print()
    except Exception as e:
        print(f"  [错误] 初始化失败: {e}")
        print("  请检查 .env 文件中的 DEEPSEEK_API_KEY 配置。")
        sys.exit(1)

    while True:
        try:
            # 多行输入
            print(f"\n{DIVIDER_THIN}")
            print("  You (多行输入, 连按两次 Enter 发送):")
            print(DIVIDER_THIN)
            lines = []
            while True:
                try:
                    line = input("  > ")
                except EOFError:
                    break
                if line.strip() == "" and lines:
                    break
                lines.append(line)

            user_input = "\n".join(lines).strip()
            if not user_input:
                continue

            # 命令处理
            if user_input.startswith("/"):
                cmd = user_input.lower().strip()
                if cmd in ("/quit", "/exit", "/q"):
                    print(f"\n  [退出] 感谢使用 Gener. 再见!\n")
                    break
                elif cmd == "/reset":
                    agent.reset()
                    continue
                elif cmd == "/tools":
                    print(f"\n{DIVIDER_THIN}")
                    print("  已注册工具列表:")
                    print(DIVIDER_THIN)
                    for i, name in enumerate(agent.get_tool_list(), 1):
                        print(f"  {i:2d}. {name}")
                    print(DIVIDER_THIN)
                    continue
                elif cmd == "/help":
                    print(HELP_TEXT)
                    continue
                else:
                    print(f"  [提示] 未知命令: {cmd}")
                    continue

            # Agent 对话（流式逐字输出）
            print(f"\n{DIVIDER_THIN}")
            print("  Gener:")
            print(DIVIDER_THIN)
            reply = agent.chat_stream(user_input)
            # reply 已在 chat_stream 内部逐字打印
            print()

        except KeyboardInterrupt:
            print(f"\n\n  [退出] 感谢使用 Gener. 再见!\n")
            break
        except Exception as e:
            print(f"\n  [错误] {e}\n")
            import traceback
            traceback.print_exc()


if __name__ == "__main__":
    main()
