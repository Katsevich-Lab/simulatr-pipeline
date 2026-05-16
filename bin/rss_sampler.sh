#!/usr/bin/env bash
# rss_sampler.sh start | stop
#
# Background RSS sampler for non-Slurm executors (laptop, Docker).
# `start` spawns a child that polls `ps -o rss= -p $PPID` once per
# second and tracks the maximum RSS in KB. `stop` kills the sampler
# and writes the result to .command.memory.tsv in a format the
# COLLECT_MEMORY process understands. Under Slurm, conf/slurm.config
# overwrites .command.memory.tsv via `sacct` in afterScript, so this
# sampler's output is harmlessly clobbered.

set -u

PARENT_PID=$$
PIDFILE=".rss_sampler.pid"
MAXFILE=".rss_sampler.max"

case "${1:-}" in
  start)
    rm -f "$MAXFILE"
    (
      max=0
      while kill -0 "$PARENT_PID" 2>/dev/null; do
        cur=$(ps -o rss= -p "$PARENT_PID" 2>/dev/null | tr -d ' ' || echo 0)
        if [ -n "$cur" ] && [ "$cur" -gt "$max" ]; then
          max=$cur
        fi
        echo "$max" > "$MAXFILE"
        sleep 1
      done
    ) &
    echo $! > "$PIDFILE"
    ;;

  stop)
    if [ -f "$PIDFILE" ]; then
      kill "$(cat "$PIDFILE")" 2>/dev/null || true
      rm -f "$PIDFILE"
    fi
    sleep 1
    max_kb=$(cat "$MAXFILE" 2>/dev/null || echo 0)
    rm -f "$MAXFILE"
    {
      echo "JobID|JobName|MaxRSS|MaxVMSize|Elapsed"
      echo "rss_sampler|local|${max_kb}K||"
    } > .command.memory.tsv
    ;;

  *)
    echo "usage: $0 start|stop" >&2
    exit 64
    ;;
esac
