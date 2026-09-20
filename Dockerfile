# 运行时锁定 Python 3.11（slim 基于 Debian bookworm）。
FROM python:3.11-slim-bookworm AS base

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /srv

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY app ./app

EXPOSE 8000

# 容器化运行 uvicorn：单进程内 async 路由即可；如需多进程可用：
#   uvicorn app.main:app --host 0.0.0.0 --port 8000 --workers 2
# 服务无共享可变状态，多 worker/多请求并发互不干扰。
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]

# ── 测试镜像目标：docker build --target test -t lame-test . && docker run --rm lame-test
FROM base AS test
COPY requirements-dev.txt ./
RUN pip install --no-cache-dir -r requirements-dev.txt
COPY conftest.py pytest.ini ./
COPY tests ./tests
CMD ["pytest", "-q"]
