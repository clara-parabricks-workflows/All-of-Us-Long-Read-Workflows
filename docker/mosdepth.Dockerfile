FROM erictdawson/base
RUN wget https://github.com/brentp/mosdepth/releases/download/v0.3.5/mosdepth && \
    mv mosdepth /usr/bin/ && \
    chmod 777 /usr/bin/mosdepth
ENTRYPOINT bash
