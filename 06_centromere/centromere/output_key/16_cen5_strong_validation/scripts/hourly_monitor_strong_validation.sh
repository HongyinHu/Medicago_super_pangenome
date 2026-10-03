#!/usr/bin/env bash
set -euo pipefail

RUN=path/to/project/10.centromere_analysis/output_key/16_cen5_strong_validation
MONITOR_LOG="${RUN}/logs/hourly_monitor.log"
PID_FILE="${RUN}/logs/nohup_strong_validation.pid"

cd "${RUN}"

status_once() {
  local now
  now=$(date '+%F %T')
  echo "===== ${now} =====" >> "${MONITOR_LOG}"

  if [[ -s "${PID_FILE}" ]]; then
    local pid
    pid=$(cat "${PID_FILE}")
    if ps -p "${pid}" >/dev/null 2>&1; then
      echo "pipeline_pid=${pid} status=RUNNING" >> "${MONITOR_LOG}"
      ps -p "${pid}" -o pid,ppid,stat,etime,cmd --no-headers >> "${MONITOR_LOG}" 2>/dev/null || true
    else
      echo "pipeline_pid=${pid} status=NOT_RUNNING" >> "${MONITOR_LOG}"
    fi
  else
    echo "pipeline_pid=NA status=NO_PID_FILE" >> "${MONITOR_LOG}"
  fi

  if grep -q 'Strong validation pipeline finished' logs/run_strong_validation.log 2>/dev/null; then
    echo "completion=FINISHED_MARKER_FOUND" >> "${MONITOR_LOG}"
  else
    echo "completion=NOT_FINISHED" >> "${MONITOR_LOG}"
  fi

  echo "latest_main_log:" >> "${MONITOR_LOG}"
  tail -20 logs/run_strong_validation.log >> "${MONITOR_LOG}" 2>/dev/null || true

  echo "recent_error_signals:" >> "${MONITOR_LOG}"
  grep -RniE 'error|failed|traceback|killed|no such|missing|cannot|segmentation' logs 2>/dev/null | tail -30 >> "${MONITOR_LOG}" || true

  echo "key_outputs:" >> "${MONITOR_LOG}"
  for f in \
    results/junction_refined_intervals.tsv \
    results/targeted_repeat_decay_summary.tsv \
    results/x8_raw_q20_called_CENH3_domains.tsv \
    results/x8_CEN5_pm3Mb_projection_to_Mpo_summary.tsv
  do
    if [[ -s "${f}" ]]; then
      printf "%s\t%s bytes\t%s\n" "${f}" "$(stat -c%s "${f}")" "$(date -r "${f}" '+%F %T')" >> "${MONITOR_LOG}"
    else
      printf "%s\tMISSING_OR_EMPTY\n" "${f}" >> "${MONITOR_LOG}"
    fi
  done

  echo "disk:" >> "${MONITOR_LOG}"
  df -h "${RUN}" >> "${MONITOR_LOG}" 2>/dev/null || true
  echo >> "${MONITOR_LOG}"
}

status_once

while true; do
  sleep 3600
  status_once
done
