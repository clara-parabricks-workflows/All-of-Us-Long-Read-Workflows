FROM nvcr.io/nvidia/clara/clara-parabricks:4.0.0-1

COPY --from=claraparabricks/bwa:latest /usr/bin/bwa /usr/bin/
COPY --from=claraparabricks/bwa:latest /usr/local/bin/bgzip /usr/bin/
COPY --from=claraparabricks/bwa:latest /usr/local/bin/tabix /usr/bin/

## Strip off the entrypoint
ENTRYPOINT bash