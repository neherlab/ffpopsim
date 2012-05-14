/* renames and ignores */
%ignore SAMPLE_ERROR;
%ignore sample;

/**** INDEX_VALUE_PAIR_T ****/
%rename(index_value_pair) index_value_pair_t;
%extend index_value_pair_t {
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"index: %d, val: %17.2e", $self->index, $self->val);
        return &buffer[0];
}

const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"(%d, %17.2e)", $self->index, $self->val);
        return &buffer[0];
}
} /* extend index_value_pair_t */

/**** STAT_T ****/
%rename(stat) stat_t;
%extend stat_t {
const char* __str__() {
        static char buffer[255];
        sprintf(buffer,"mean: %17.2e, variance: %17.2e", $self->mean, $self->variance);
        return &buffer[0];
}

const char* __repr__() {
        static char buffer[255];
        sprintf(buffer,"(%17.2e, %17.2e)", $self->mean, $self->variance);
        return &buffer[0];
}
} /* extend stat_t */
