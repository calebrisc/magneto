/* M1 bench logger.
 *
 * Streams two differential SAADC channels (sin/cos TMR bridge pair)
 * as CSV over the UART console:
 *
 *   t_ms,ch0_raw,ch0_mv,ch1_raw,ch1_mv
 *
 * SAMPLE_HZ is capped by console bandwidth: ~30 bytes/line at
 * 115200 baud tops out near 380 lines/s. 200 Hz default leaves
 * headroom; gap characterization is static/slow-slide data anyway.
 * For flick-speed captures later, raise the VCOM baud and this rate.
 */

#include <zephyr/kernel.h>
#include <zephyr/drivers/adc.h>
#include <zephyr/sys/printk.h>

#define SAMPLE_HZ 200

static const struct adc_dt_spec channels[] = {
	ADC_DT_SPEC_GET_BY_IDX(DT_PATH(zephyr_user), 0),
	ADC_DT_SPEC_GET_BY_IDX(DT_PATH(zephyr_user), 1),
};

int main(void)
{
	int16_t buf;
	struct adc_sequence sequence = {
		.buffer = &buf,
		.buffer_size = sizeof(buf),
	};

	for (size_t i = 0; i < ARRAY_SIZE(channels); i++) {
		if (!adc_is_ready_dt(&channels[i])) {
			printk("# ADC %s not ready\n", channels[i].dev->name);
			return 0;
		}
		int err = adc_channel_setup_dt(&channels[i]);
		if (err) {
			printk("# channel %d setup failed (%d)\n", (int)i, err);
			return 0;
		}
	}

	printk("# m1_logger: %d Hz, differential, 12-bit\n", SAMPLE_HZ);
	printk("t_ms,ch0_raw,ch0_mv,ch1_raw,ch1_mv\n");

	while (1) {
		int64_t t = k_uptime_get();
		int32_t raw[ARRAY_SIZE(channels)];
		int32_t mv[ARRAY_SIZE(channels)];

		for (size_t i = 0; i < ARRAY_SIZE(channels); i++) {
			adc_sequence_init_dt(&channels[i], &sequence);
			int err = adc_read_dt(&channels[i], &sequence);
			if (err) {
				printk("# read ch%d failed (%d)\n", (int)i, err);
				raw[i] = 0;
				mv[i] = 0;
				continue;
			}
			raw[i] = buf;
			mv[i] = buf;
			adc_raw_to_millivolts_dt(&channels[i], &mv[i]);
		}

		printk("%lld,%d,%d,%d,%d\n", t, raw[0], mv[0], raw[1], mv[1]);
		k_sleep(K_USEC(1000000 / SAMPLE_HZ));
	}
	return 0;
}
